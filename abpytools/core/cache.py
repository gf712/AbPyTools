class Cache:
    def __init__(self, max_cache_size=10):
        self.max_cache_size = max_cache_size
        self.cache = {}
        self._cache_size = 0
        self._cache_keys = []

    def update(self, key, data, override=True):
        # make space if the cache is full and override is True
        if len(self) >= self.max_cache_size and override:
            self._clear_cache()

        if len(self) >= self.max_cache_size and not override:
            raise ValueError("Cache is full, either increase the size of cache, "
                             "remove items from cache or allow override")

        if key not in self:
            self.add(key, data)

    def _clear_cache(self):
        while len(self) >= self.max_cache_size:
            # clear up cache until reaching max_cache_size - 1
            self.remove(self._cache_keys[0])

    def remove(self, key):
        self.cache.pop(key)
        self._cache_keys.pop(self._cache_keys.index(key))
        self._cache_size -= 1

    def add(self, key, data):
        self.cache[key] = data
        self._cache_keys.append(key)
        self._cache_size += 1

    def empty_cache(self):
        self.cache = {}
        self._cache_keys = []
        self._cache_size = 0

    def _string_summary_basic(self):
        return "abpytools.Cache size: {}".format(len(self))

    def __getitem__(self, item):
        return self.cache[item]

    def __contains__(self, item):
        if item in self._cache_keys:
            return True
        else:
            return False

    def __len__(self):
        return self._cache_size

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))

    def __setitem__(self, key, value):
        self.update(key, value)
