import unittest

def print_suite(suite):
    tests = []

    if hasattr(suite, '__iter__'):
        for x in suite:
            print_suite(x)
    else:
        tests.append(suite.id())

    tests.sort()
    for t in tests:
        print(t)

if __name__ == '__main__':
    print_suite(unittest.defaultTestLoader.discover('.'))
