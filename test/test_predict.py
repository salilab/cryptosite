import unittest
import utils
import os
import sys
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
    std_ = None

    def transform(self, inp):
        return numpy.array([[-1.48500907, -0.98061114, -0.3578966],
                            [-0.4486418, -0.90488697, -0.1571045]])


class MockSVM(object):
    impl = None
    class_weight_label_ = None
    dual_coef_ = None

    def predict(self, inp):
        return numpy.array([0., 0.])

    def predict_proba(self, inp):
        return numpy.array([[0.97945919, 0.02054081],
                            [0.98243048, 0.01756952]])


def mock_pickle_loads(data, encoding=None):
    if b'Scaler' in data:
        return MockScaler()
    else:
        return MockSVM()


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to predict"""
        for args in (['x', 'y'], []):
            _ = utils.check_output(['cryptosite', 'predict'] + args,
                                   stderr=subprocess.STDOUT, retcode=2)
            _ = utils.check_output([sys.executable, '-m',
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

    def test_main_mocked(self):
        """Test complete run of predict (with mocked sklearn)"""
        indir = os.path.join(TOPDIR, 'test', 'input')
        with utils.temporary_working_directory():
            shutil.copy(os.path.join(indir, 'XXX.features'), '.')
            shutil.copy(os.path.join(indir, 'XXX_mdl.pdb'), '.')
            with utils.mocked_object(pickle, 'loads', mock_pickle_loads):
                cryptosite.predict.predict('XXX.features', model='final')
            os.unlink('XXX.pol.pred')
            os.unlink('XXX.pol.pred.pdb')

    def test_main_real(self):
        """Test complete run of predict (with real sklearn)"""
        indir = os.path.join(TOPDIR, 'test', 'input')
        with utils.temporary_working_directory():
            shutil.copy(os.path.join(indir, 'XXX.features'), '.')
            shutil.copy(os.path.join(indir, 'XXX_mdl.pdb'), '.')
            cryptosite.predict.predict('XXX.features', model='final')
            with open('XXX.pol.pred') as fh:
                lines = fh.readlines()
            lines = sorted(line for line in lines
                           if not line.startswith('PDBID'))
            scores = [line.rstrip('\r\n').split('\t') for line in lines]
            os.unlink('XXX.pol.pred')
            os.unlink('XXX.pol.pred.pdb')
        # Compare scores against known-good values
        sqc, ptm, cnc, crypt = [float(x) for x in scores[0][3:]]
        self.assertAlmostEqual(sqc, -1.4850, delta=1e-4)
        self.assertAlmostEqual(ptm, -0.9806, delta=1e-4)
        self.assertAlmostEqual(cnc, -0.3579, delta=1e-4)
        self.assertAlmostEqual(crypt, 0.0205, delta=1e-4)

        sqc, ptm, cnc, crypt = [float(x) for x in scores[1][3:]]
        self.assertAlmostEqual(sqc, -0.4486, delta=1e-4)
        self.assertAlmostEqual(ptm, -0.9049, delta=1e-4)
        self.assertAlmostEqual(cnc, -0.1571, delta=1e-4)
        self.assertAlmostEqual(crypt, 0.0176, delta=1e-4)


if __name__ == '__main__':
    unittest.main()
