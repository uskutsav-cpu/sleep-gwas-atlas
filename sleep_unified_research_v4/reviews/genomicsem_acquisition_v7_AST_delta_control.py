#!/usr/bin/env python3
"""Exact AST allowlist for v7 origin import, receipt path and count changes."""
import ast
import copy
import hashlib
import json
from pathlib import Path

P=Path(__file__).resolve().parents[1]
OLD=P/'scripts/52_acquire_extension_raw_sources_v6.py'
NEW=P/'scripts/52_acquire_extension_raw_sources_v7.py'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
expected={str(OLD):'d408b245c362c274e95f8c33727ba016367ed14b1e0cd0a77adfd29ae577bf8d',str(NEW):'f4f5e7b8c5c1f4942ff65a95737b525c1df9f9f53633c1fa269a1c774976de12'}
for p,h in expected.items():assert sha(p)==h,p
functions=lambda p:{n.name:n for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
a,b=functions(OLD),functions(NEW)

class AllowedDelta(ast.NodeTransformer):
    def visit_Constant(self,node):
        if isinstance(node.value,str):
            node.value=node.value.replace('extension_replay_common_v3.py','extension_replay_common_v2.py').replace('extension_raw_acquisition_family_receipt_v4_7.json','extension_raw_acquisition_family_receipt_v4_6.json')
        return node
    def visit_Call(self,node):
        node=self.generic_visit(node)
        if isinstance(node.func,ast.Name) and node.func.id=='len' and len(node.args)==1 and ast.unparse(node.args[0])=="plan['reused_checkpoints']":
            return ast.copy_location(ast.Constant(value=2),node)
        return node
    def visit_BinOp(self,node):
        node=self.generic_visit(node)
        if isinstance(node.op,ast.Sub) and isinstance(node.left,ast.Constant) and node.left.value==100 and isinstance(node.right,ast.Constant) and node.right.value==2:
            return ast.copy_location(ast.Constant(value=98),node)
        return node

checks={}
for name in ['reused_checkpoint_gate','finalize_family','execute']:
    checks[name]=ast.dump(a[name],include_attributes=False)==ast.dump(AllowedDelta().visit(copy.deepcopy(b[name])),include_attributes=False)
    assert checks[name],name
for p,h in expected.items():assert sha(p)==h,p
out=P/'reviews/genomicsem_acquisition_v7_AST_delta_receipt.json'
with out.open('x') as f:f.write(json.dumps({'status':'EXACT_AST_DELTA_ALLOWLIST_PASS','source_sha256_before_and_after':expected,'checks':checks,'normalizations_only':['common origin gate v3 import mapped to v2','family receipt version7 path mapped to6','len(plan.reused_checkpoints) mapped to prior exact2;100-minus-length mapped to98'],'helper_sha256':sha(__file__),'body_reads':0,'workers':0},indent=2)+'\n')
print(json.dumps(checks))
