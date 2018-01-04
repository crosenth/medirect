import unittest

from medirect import ftract


class Test_ftract(unittest.TestCase):

    def setUp(self):
        with open('tests/data.ft') as data:
            self.data = list(data)

    def test01(self):
        correct = [('1', 223523, 225078, '2'), ('629', 150, 1687, '1')]
        results = ftract.Ftract.filter_features(
            None, self.data, ['rrna:product:16S'])
        self.assertEqual(correct, list(results))

    def test02(self):
        correct = [
            ('1', 304, 1728, '1'),
            ('1', 1866, 2657, '1'),
            ('1', 2650, 3582, '1'),
            ('1', 3560, 3769, '1'),
            ('1', 3773, 4867, '1'),
            ('1', 4848, 6149, '1'),
            ('1', 6152, 6412, '1'),
            ('1', 223523, 225078, '2'),
            ('629', 150, 1687, '1')]
        results = ftract.Ftract.filter_features(
            None, self.data, ['::'])
        self.assertEqual(correct, list(results))

    def test03(self):
        correct = [('1', 223523, 225078, '2'), ('629', 150, 1687, '1')]
        results = ftract.Ftract.filter_features(
            None, self.data, ['rrna::'])
        self.assertEqual(correct, list(results))

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
            ('629', 150, 1687, '1')]
        results = ftract.Ftract.filter_features(
            None, self.data, [':product:'])
        self.assertEqual(correct, list(results))

    def test05(self):
        correct = [
            ('1', 304, 1728, '1'),
            ('1', 1866, 2657, '1'),
            ('1', 2650, 3582, '1'),
            ('1', 3560, 3769, '1'),
            ('1', 3773, 4867, '1'),
            ('1', 4848, 6149, '1'),
            ('1', 6152, 6412, '1'),
            ('1', 223523, 225078, '2'),
            ('629', 150, 1687, '1')]
        results = ftract.Ftract.filter_features(
            None, self.data, [':product:'])
        self.assertEqual(correct, list(results))

    def test06(self):
        correct = [('1', 3560, 3769, '1')]
        results = ftract.Ftract.filter_features(
            None, self.data, ['::OJPFPCPC_00004'])
        self.assertEqual(correct, list(results))


if __name__ == '__main__':
    unittest.main()
