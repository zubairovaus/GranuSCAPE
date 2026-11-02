"""Monotonic elapsed-time measurements; nested spans are not additive."""

from contextlib import contextmanager
from time import perf_counter


class Timer:
    def __init__(self):
        self.seconds = {}

    @contextmanager
    def span(self, name):
        start = perf_counter()
        try:
            yield
        finally:
            self.seconds[name] = self.seconds.get(name, 0.0) + perf_counter() - start
