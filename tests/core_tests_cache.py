import unittest
from abpytools.core.cache import Cache


class CacheTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = [[0, 1], 'n', ('tuple', 0), 3.1415]
        cls.names = ['a', 'b', 'c', 'd']
        cls.data_dict = {key: value for key, value in zip(cls.names, cls.data)}

    def test_cache_instantiation(self):
        cache = Cache()
        self.assertIsInstance(cache, Cache)

    def test_cache_update_1(self):
        cache = Cache()
        cache.update(key='a', data=self.data_dict['a'])
        self.assertEqual(cache['a'], self.data_dict['a'])

    def test_cache_update_2(self):
        # check if cache works when adding a second item
        cache = Cache()
        cache.update(key='a', data=self.data_dict['a'])
        cache.update(key='b', data=self.data_dict['b'])
        self.assertEqual(cache['b'], self.data_dict['b'])

    def test_contains(self):
        cache = Cache()
        cache.update(key='a', data=self.data_dict['a'])
        self.assertTrue('a' in cache)

    def test_cache_remove_key(self):
        cache = Cache()
        cache.update(key='a', data=self.data_dict['a'])
        cache.update(key='b', data=self.data_dict['b'])
        cache.remove('a')
        self.assertFalse('a' in cache)

    def test_cache_empty(self):
        cache = Cache()
        cache.update(key='a', data=self.data_dict['a'])
        cache.update(key='b', data=self.data_dict['b'])
        cache.empty_cache()
        self.assertEqual(len(cache), 0)

    def test_cache_maximum_size_1(self):
        cache = Cache(max_cache_size=2)
        cache.update(key='a', data=self.data_dict['a'])
        cache.update(key='b', data=self.data_dict['b'])
        # override is True by default, but this makes it explicit that it is testing for this (default) situation
        cache.update(key='c', data=self.data_dict['c'], override=True)
        self.assertFalse('a' in cache)

    def test_cache_maximum_size_2(self):
        cache = Cache(max_cache_size=2)
        cache.update(key='a', data=self.data_dict['a'])
        cache.update(key='b', data=self.data_dict['b'])
        cache.update(key='c', data=self.data_dict['c'], override=True)
        self.assertEqual(len(cache), 2)

    def test_cache_maximum_size_exception(self):
        cache = Cache(max_cache_size=2)
        cache.update(key='a', data=self.data_dict['a'])
        cache.update(key='b', data=self.data_dict['b'])
        self.assertRaises(ValueError, cache.update, 'c', self.data_dict['c'], False)
