from collections import Counter
from heapq import heappop, heappush


class BoundedItemsHeap:
    # Collection of items with lower/upper bounds on their "values" (item.lo, item.up)
    #
    # Two thresholds may be set:
    # * good - if item.up <= good then item is considered "good" and discarded
    # * bad - if item.lo >= bad then item is "bad" and it is temporarily stored
    # note that an item may be good and bad simultaneously, if bad_thr <= lo <= up <= good_thr
    # in this case it is considered as "good"
    # * otherwise item is "active", they are stored in heap with priority=up

    def __init__(self, keep_max_lo_item=False, stash=None):
        # use heapq algorithm for list of tuples (priority, increment, item)
        # keep_max_lo_item: subj, note that item may become non-active
        # stash - to store some per-instance data
        self._heap = []
        self._bad_items = []
        self._good_threshold = None
        self._bad_threshold = None

        self._keep_max_lo_item = keep_max_lo_item
        if self._keep_max_lo_item:
            self.max_lo_item = None

        self.stats = Counter()
        self.stash = stash
        self._inc = 0

    def set_good_threshold(self, threshold):
        # If up <= good, item is considered "good" (hence insignificant).
        self.stats['set_good_threshold'] += 1
        self._good_threshold = threshold

    def set_bad_threshold(self, threshold):
        # If lo >= bad, item is considered "bad".
        self.stats['set_bad_threshold'] += 1
        self._bad_threshold = threshold

    def has_items(self):
        return bool(self._heap)

    def active_items(self):
        # returns iterator over active items
        return (node[-1] for node in self._heap)

    def push(self, item):
        # Add node checking the thresholds:
        # good pair is dropped / bad pair is temporarily stored / otherwise node is added to the heap
        self.stats['push'] += 1

        if self._keep_max_lo_item:
            if self.max_lo_item is None or item.lo > self.max_lo_item.lo:
                self.max_lo_item = item

        # the order of checks may be important
        # is Estimator we add bad pairs to SAT clauses, so it may be more effective to
        # check for goodness first (recall that an item may be good and bad simultaneously)
        if (self._good_threshold is not None) and item.up <= self._good_threshold:
            self.stats['good'] += 1
            return
        if (self._bad_threshold is not None) and item.lo >= self._bad_threshold:
            self._bad_items.append(item)
            self.stats['bad'] += 1
            return

        self._inc += 1
        node = (-item.up, self._inc, item)  # first key is priority; _inc to avoid comparing items
        heappush(self._heap, node)

    def _extend(self, items):
        for item in items:
            self.push(item)

    def top(self):
        # Active item with highest priority (up)
        return self._heap[0][-1]

    def pop(self):
        # Pop and return active item with highest priority (up)
        return heappop(self._heap)[-1]

    def pop_bad_items(self):
        # Return bad pairs list and empty it.
        bad_items = self._bad_items
        self._bad_items = []
        return bad_items

    def cleanup(self):
        # check that cleanup pushs <= regular pushs to avoid too much cleanups
        if self.stats['cleanup_push'] <= (self.stats['push'] - self.stats['cleanup_push']):
            self._do_cleanup()

    def _do_cleanup(self):
        # rebuild heap with actual thresholds
        items = list(self.active_items())
        self._heap = []
        self._extend(items)
        self.stats['cleanup_push'] += len(items)
        self.stats['cleanup_count'] += 1

    def copy_and_cleanup(self, good_threshold=None, bad_threshold=None):
        # copy initial tree and apply thresholds, if any
        assert self._good_threshold is None
        assert self._bad_threshold is None
        assert not self._keep_max_lo_item
        new_heap = BoundedItemsHeap()
        new_heap.set_good_threshold(good_threshold)
        new_heap.set_bad_threshold(bad_threshold)
        old_items = list(self.active_items())
        new_heap._extend(old_items)
        new_heap.stats['copy_push'] += len(old_items)
        return new_heap
