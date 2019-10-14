import unittest
import utils
import os
import sys
import re
import subprocess
import shutil
import pickle
import numpy

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)
import cryptosite.predict

# Replace the pickled scaler and SVM objects with mocks. These will run
# faster, and don't require that our test machine has the exact same version
# of sklearn as the production machine (pickles aren't portable to different
# versions of sklearn).
class MockScaler(object):
    def transform(self, inp):
        return numpy.array([[-1.48500907, -0.98061114, -0.3578966 ],
                            [-0.4486418 , -0.90488697, -0.1571045 ]])

class MockSVM(object):
    def predict(self, inp):
        return numpy.array([ 0.,  0.])

    def predict_proba(self, inp):
        return numpy.array([[ 0.97945919,  0.02054081],
                            [ 0.98243048,  0.01756952]])

def mock_pickle_load(fh):
    if 'Scaler' in fh.read():
        return MockScaler()
    else:
        return MockSVM()

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to predict"""
        for args in (['x', 'y'], []):
            out = utils.check_output(['cryptosite', 'predict'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)
            out = utils.check_output([sys.executable, '-m',
                                     'cryptosite.predict'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)

    def test_get_matrix(self):
        """Test get_matrix() function"""
        fname = os.path.join(TOPDIR, 'test', 'input', 'test.features')
        m, header, indices = cryptosite.predict.get_matrix(fname,
                                                           model='final')
        self.assertEqual(len(m), 2)
        self.assertEqual(header, ['SQC', 'PTM', 'CNC_mean_'])
        self.assertEqual(len(indices), 2)

    def test_main(self):
        """Test complete run of predict"""
        indir = os.path.join(TOPDIR, 'test', 'input')
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(indir, 'XXX.features'), '.')
            shutil.copy(os.path.join(indir, 'XXX_mdl.pdb'), '.')
            with utils.mocked_object(pickle, 'load', mock_pickle_load):
                cryptosite.predict.predict('XXX.features', model='final')
            os.unlink('XXX.pol.pred')
            os.unlink('XXX.pol.pred.pdb')

if __name__ == '__main__':
    unittest.main()
