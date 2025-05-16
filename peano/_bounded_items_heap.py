from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from heapq import heappop, heappush
from typing import Any, Iterable
from numbers import Number
from quicktions import Fraction  # type: ignore


@dataclass
class _HeapStat:
    update_good_threshold: int = 0
    rebuild_count: int = 0
    push: int = 0
    copy_push: int = 0
    good: int = 0
    bad: int = 0


@dataclass(frozen=True)
class BoundedItem:
    lo: Fraction
    up: Fraction


class BoundedItemsHeap[T: BoundedItem]:
    # Collection of items with lower/upper bounds on their "values" (item.lo, item.up)
    # T = items type
    #
    # Two thresholds may be set:
    # * good - if item.up <= good then item is considered "good" and discarded
    # * bad - if item.lo >= bad then item is "bad" and it is temporarily stored
    # note that an item may be good and bad simultaneously, if bad_thr <= lo <= up <= good_thr
    # in this case it is considered as "good"
    # * otherwise item is "active", they are stored in heap with priority=up

    # TODO: get rid of stash

    _heap: list[tuple[Fraction, int, T]]
    _bad_items: list[T]
    max_lo_item: T

    def __init__(self, good_threshold: Fraction | None, bad_threshold: Fraction | None, keep_max_lo_item: bool = False, stash: Any = None) -> None:
        # use heapq algorithm for nodes in _heap (priority, increment, item)
        # keep_max_lo_item: subj, note that item may become non-active
        # stash - to store some per-instance data
        self._heap = []
        self._bad_items = []
        self._good_threshold = good_threshold
        self._bad_threshold = bad_threshold

        self._keep_max_lo_item = keep_max_lo_item
        if self._keep_max_lo_item:
            self.seen_items = False

        self.stats = _HeapStat()
        self.stash = stash
        self._inc = 0

    def update_good_threshold(self, threshold: Fraction) -> None:
        # we can only raise threshold to keep heap invariant correct
        if (self._good_threshold is not None) and threshold < self._good_threshold:
            raise ValueError("Can't set lower good threshold!")

        self._good_threshold = threshold
        self.stats.update_good_threshold += 1

        # cleanup to maintain heap invariant
        active_items = [node[-1] for node in self._heap if node[-1].up > threshold]
        if len(active_items) < len(self._heap):
            self._heap = []
            self._extend(active_items)
            self.stats.rebuild_count += 1

    def size(self) -> int:
        return len(self._heap)

    def items(self) -> Iterable[T]:
        for node in self._heap:
            yield node[-1]

    def push(self, item: T) -> None:
        # Add node checking the thresholds:
        # good pair is dropped / bad pair is temporarily stored / otherwise node is added to the heap
        self.stats.push += 1

        if self._keep_max_lo_item:
            if (not self.seen_items) or item.lo > self.max_lo_item.lo:
                self.max_lo_item = item
                self.seen_items = True

        # the order of checks may be important
        # is Estimator we add bad pairs to SAT clauses, so it may be more effective to
        # check for goodness first (recall that an item may be good and bad simultaneously)
        if (self._good_threshold is not None) and item.up <= self._good_threshold:
            self.stats.good += 1
            return
        if (self._bad_threshold is not None) and item.lo >= self._bad_threshold:
            self._bad_items.append(item)
            self.stats.bad += 1
            return

        self._inc += 1
        node = (-item.up, self._inc, item)  # first key is priority; _inc to avoid comparing items
        heappush(self._heap, node)

    def _extend(self, items: Iterable[T]) -> None:
        # one may use append+heapify instead of heappush, but this is not faster
        for item in items:
            self.push(item)

    def top(self) -> T:
        # Active item with highest priority (up)
        return self._heap[0][-1]

    def pop(self) -> T:
        # Pop and return active item with highest priority (up)
        return heappop(self._heap)[-1]

    def pop_bad_items(self) -> list[T]:
        # Return bad pairs list and empty it.
        bad_items = self._bad_items
        self._bad_items = []
        return bad_items

    def copy_and_cleanup(self, good_threshold: Fraction | None = None, bad_threshold: Fraction | None = None) -> BoundedItemsHeap[T]:
        # copy initial tree and apply thresholds, if any
        assert self._good_threshold is None
        assert self._bad_threshold is None
        assert not self._keep_max_lo_item
        new_heap: BoundedItemsHeap[T] = BoundedItemsHeap(good_threshold, bad_threshold)
        new_heap._extend(self.items())
        new_heap.stats.copy_push += self.size()
        return new_heap
