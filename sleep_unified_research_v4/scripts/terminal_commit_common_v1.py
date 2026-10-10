"""Durable, fail-closed metadata commit for a reviewed research stage.

All result receipts are provisional while PENDING exists. A caller must hold
its existing ownership/cleanup boundary and supply the independently reviewed
identity, resource and termination gates. This module never launches a worker.
"""
import hashlib
import json
import os
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def save_new(path, value):
    payload = json.dumps(value, indent=2, allow_nan=False) + '\n'
    with Path(path).open('x') as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())


def sync_directory(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


class TerminalCommit:
    def __init__(self, pending, seal, binding):
        self.pending, self.seal = Path(pending), Path(seal)
        self.binding = dict(binding)
        if self.pending.parent != self.seal.parent or self.pending.parent.is_symlink():
            raise RuntimeError('ONE_REGULAR_PRIVATE_TERMINAL_DIRECTORY_REQUIRED')
        if self.pending.exists() or self.pending.is_symlink() or self.seal.exists() or self.seal.is_symlink():
            raise RuntimeError('PRIOR_TERMINAL_ATTEMPT_PRESERVED_NO_RETRY')
        save_new(self.pending, dict(status='PENDING_NOT_ADMISSIBLE', binding=self.binding, owner_pid=os.getpid()))
        sync_directory(self.pending.parent)
        self.pending_sha = sha(self.pending)

    def commit(self, receipt_sha256, *, identity_gate, resource_gate, termination_gate):
        """Atomic final marker removal follows all persistence and checks.

        A failed commit leaves PENDING, even if supplemental failure recording
        also fails. Nothing fallible follows successful marker removal.
        """
        expected = dict(receipt_sha256)
        try:
            if not expected or any(not Path(p).is_file() or Path(p).is_symlink() or sha(p) != h for p, h in expected.items()):
                raise RuntimeError('EXACT_NONEMPTY_RESULT_RECEIPT_BINDINGS_REQUIRED')
            if not self.pending.is_file() or sha(self.pending) != self.pending_sha:
                raise RuntimeError('DURABLE_PENDING_CHANGED_OR_MISSING')
            identity_gate()
            termination_gate()
            resource_gate()
            save_new(self.seal, dict(status='REVIEWED_STAGE_TERMINAL_SEAL', binding=self.binding,
                                    result_receipt_sha256=expected, pending_path=str(self.pending),
                                    pending_sha256=self.pending_sha,
                                    success_requires_absent_PENDING_and_all_failure_addenda=True))
            sync_directory(self.seal.parent)
            seal_sha = sha(self.seal)
            identity_gate()
            if any(sha(p) != h for p, h in expected.items()) or sha(self.seal) != seal_sha or sha(self.pending) != self.pending_sha:
                raise RuntimeError('TERMINAL_RESULT_OR_PENDING_CHANGED_AFTER_PERSISTENCE')
            for p in [*expected, str(self.seal)]:
                if Path(p + '.failure.json').exists():
                    raise RuntimeError('TERMINAL_FAILURE_INVALIDATES_PROVISIONAL_RECEIPTS')
            resource_gate()
            termination_gate()
            self.pending.unlink()
            return True
        except BaseException as error:
            try:
                save_new(Path(str(self.seal) + '.failure.json'),
                         dict(status='TERMINAL_COMMIT_FAILED_PRESERVED', binding=self.binding,
                              error=type(error).__name__ + ': ' + str(error), PENDING_remains_authoritative_veto=True))
            except BaseException:
                pass
            return False


def require_committed(pending, seal, binding, receipt_sha256):
    """Read-only consumer; no publication credit from a provisional receipt."""
    pending, seal = Path(pending), Path(seal)
    if pending.exists() or pending.is_symlink() or not seal.is_file() or seal.is_symlink():
        raise RuntimeError('STAGE_TERMINAL_COMMIT_NOT_COMPLETE')
    seal_sha = sha(seal)
    s = json.loads(seal.read_text())
    expected = dict(receipt_sha256)
    if s.get('status') != 'REVIEWED_STAGE_TERMINAL_SEAL' or s.get('binding') != binding or s.get('result_receipt_sha256') != expected or s.get('pending_path') != str(pending) or not expected:
        raise RuntimeError('EXACT_STAGE_TERMINAL_SEAL_BINDINGS_REQUIRED')
    for p, h in expected.items():
        if Path(p + '.failure.json').exists() or not Path(p).is_file() or Path(p).is_symlink() or sha(p) != h:
            raise RuntimeError('STAGE_RECEIPT_CHANGED_OR_INVALIDATED')
    if Path(str(seal) + '.failure.json').exists() or sha(seal) != seal_sha or pending.exists():
        raise RuntimeError('STAGE_COMMIT_CHANGED_DURING_CONSUMPTION')
    return seal_sha
