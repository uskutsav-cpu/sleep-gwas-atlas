# Local SciPy sparse-array compatibility for Scanpy's getnnz call.
import scipy.sparse as sp, runpy, sys
if not hasattr(sp.csr_array, 'getnnz'):
 sp.csr_array.getnnz = lambda self, axis=None: self._getnnz(axis=axis)
test=sp.csr_array([[0,1],[2,3]])
assert list(test.getnnz(axis=1))==list(sp.csr_matrix(test).getnnz(axis=1))
script='work/stage4_env/bin/scdrs';sys.argv[0]=script
runpy.run_path(script,run_name='__main__')
