"""Real owned-process and guard checks; these are not biological validation."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import time
import unittest

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/30_prepare_and_run_ssd_native_campaign.py'
spec=importlib.util.spec_from_file_location('native_monitor_v4',SCRIPT)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

class MonitorChecks(unittest.TestCase):
    def test_final_internal_floor_is_three_gib(self):
        state={'internal_free_bytes':3*1024**3,'ssd_free_bytes':5*1024**3}
        self.assertIsNone(mod.final_limits(state,0,0,0))
        state['internal_free_bytes']-=1
        self.assertEqual(mod.final_limits(state,0,0,0),'INTERNAL_FULL_NATIVE_FLOOR_REACHED')

    def test_final_output_and_deadline_checks(self):
        state={'internal_free_bytes':4*1024**3,'ssd_free_bytes':6*1024**3}
        self.assertEqual(mod.final_limits(state,0,mod.OUTPUT_LIMIT+1,0),'NATIVE_OUTPUT_LIMIT')
        self.assertEqual(mod.final_limits(state,0,0,mod.STAGE_SECONDS+1),'STAGE_DEADLINE')
        self.assertEqual(mod.final_limits(state,mod.RSS_LIMIT+1,0,0),'WORKER_RSS_LIMIT')

    def test_orphaned_signal_ignoring_descendant_is_removed(self):
        code='import os,signal,time; p=os.fork(); signal.signal(signal.SIGTERM,signal.SIG_IGN) if p==0 else None; time.sleep(120) if p==0 else time.sleep(.2)'
        proc=subprocess.Popen([sys.executable,'-B','-c',code],start_new_session=True,
                              stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        try:
            proc.wait(timeout=5)
            self.assertTrue(mod.group_members(proc.pid))
            result=mod.terminate_owned(proc)
            self.assertIn('SIGKILL',result['signals'])
            self.assertEqual(result['remaining_group_members'],[])
        finally:
            if mod.group_members(proc.pid):mod.terminate_owned(proc)

if __name__=='__main__':unittest.main()
