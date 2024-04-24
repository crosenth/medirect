import ftract
import unittest


class Test_ftract(unittest.TestCase):

    defaults = dict(on_error='continue', full_format=False, min_length=None)

    def setUp(self):
        self.ftract = ftract.Ftract('test')
        with open('tests/data.ft') as data:
            self.data = list(data)

    def test01(self):
        correct = [('1', 223523, 225078, '2'), ('629', 150, 1687, '1')]
        results = self.ftract.filter_features(
            self.data, ['rrna:product:16S'], **self.defaults)
        self.assertEqual(correct, list(results))

    def test02(self):
        correct = set([
            ('1', 304, 1728, '1'),
            ('1', 1866, 2657, '1'),
            ('1', 2650, 3582, '1'),
            ('1', 3560, 3769, '1'),
            ('1', 3773, 4867, '1'),
            ('1', 4848, 6149, '1'),
            ('1', 6152, 6412, '1'),
            ('1', 223523, 225078, '2'),
            ('629', 150, 1687, '1'),
            ('gb|PKKU01000069.1|', 82536, 83733, '2'),
            ('gb|PKKU01000069.1|', 83450, 83733, '2'),
            ('gb|PKKU01000069.1|', 82536, 83451, '2'),
            ('gb|PKKU01000069.1|', 83959, 84030, '1')])
        results = self.ftract.filter_features(
            self.data, ['::'], **self.defaults)
        self.assertEqual(correct, set(results))

    def test03(self):
        correct = set([('1', 223523, 225078, '2'), ('629', 150, 1687, '1')])
        results = self.ftract.filter_features(
            self.data, ['rrna::'], **self.defaults)
        self.assertEqual(correct, set(results))

    def test04(self):
        correct = [
            ('1', 304, 1728, '1'),
            ('1', 1866, 2657, '1'),
            ('1', 2650, 3582, '1'),
            ('1', 3560, 3769, '1'),
            ('1', 3773, 4867, '1'),
            ('1', 4848, 6149, '1'),
            ('1', 6152, 6412, '1'),
            ('1', 223523, 225078, '2'),
            ('629', 150, 1687, '1'),
            ('gb|PKKU01000069.1|', 82536, 83451, '2'),
            ('gb|PKKU01000069.1|', 83959, 84030, '1')]
        results = self.ftract.filter_features(
            self.data, [':product:'], **self.defaults)
        self.assertEqual(correct, list(results))

    def test05(self):
        correct = set([
            ('1', 304, 1728, '1'),
            ('1', 1866, 2657, '1'),
            ('1', 2650, 3582, '1'),
            ('1', 3560, 3769, '1'),
            ('1', 3773, 4867, '1'),
            ('1', 4848, 6149, '1'),
            ('1', 6152, 6412, '1'),
            ('1', 223523, 225078, '2'),
            ('629', 150, 1687, '1'),
            ('gb|PKKU01000069.1|', 82536, 83733, '2'),
            ('gb|PKKU01000069.1|', 83450, 83733, '2'),
            ('gb|PKKU01000069.1|', 82536, 83451, '2'),
            ('gb|PKKU01000069.1|', 83959, 84030, '1')])
        results = self.ftract.filter_features(self.data, None, **self.defaults)
        self.assertEqual(correct, set(results))

    def test06(self):
        correct = [('1', 3560, 3769, '1')]
        results = self.ftract.filter_features(
            self.data, ['::OJPFPCPC_00004'], **self.defaults)
        self.assertEqual(correct, list(results))

    def test07(self):
        self.assertEqual((5, 8, '1'), self.ftract.coordinates(5, 8, '1'))

    def test08(self):
        self.assertEqual((5, 8, '2'), self.ftract.coordinates(8, 5, '1'))

    def test09(self):
        with self.assertRaises(ValueError):
            list(self.ftract.filter_features(
                ['invalid'], None, 'halt', False, None))


if __name__ == '__main__':
    unittest.main()
