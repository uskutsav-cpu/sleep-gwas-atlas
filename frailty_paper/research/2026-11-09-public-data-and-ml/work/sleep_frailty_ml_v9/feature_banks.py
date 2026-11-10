"""Reuse frozen public-model features without copying the large V8 banks."""
from v9_common import *


class SplitBank:
    def __init__(self, original, extra):
        self.old = np.load(original, mmap_mode='r'); self.extra = np.load(extra, mmap_mode='r')
        self.offset = len(self.old); self.shape = (self.offset+len(self.extra),)+self.old.shape[1:]
        assert self.old.shape[1:] == self.extra.shape[1:]

    def __getitem__(self, ids):
        ids = np.asarray(ids); scalar = ids.ndim == 0; flat = ids.reshape(-1)
        assert np.all((flat >= 0)&(flat < self.shape[0]))
        out = np.empty((len(flat),)+self.old.shape[1:], dtype=self.old.dtype)
        old = flat < self.offset
        out[old] = self.old[flat[old]]; out[~old] = self.extra[flat[~old]-self.offset]
        return out[0] if scalar else out.reshape(ids.shape+self.old.shape[1:])


def global_bank(): return SplitBank(V8/'pretrained_pooled128.npy', DATA/'extra_global128.npy')
def local_bank(): return SplitBank(V8/'pretrained_variant_local128.npy', DATA/'extra_local128.npy')
def reference_tokens(): return SplitBank(V8/'tokens_reference150.npy', DATA/'extra_tokens_reference150.npy')
def alternate_tokens(): return SplitBank(V8/'tokens_alternate150.npy', DATA/'extra_tokens_alternate150.npy')
def avi_bank(): return SplitBank(V8/'avi_training_bank.npy', DATA/'extra_AVI.npy')
