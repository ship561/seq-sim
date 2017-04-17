from itertools import islice


def distinct(coll):
    ''' Originally proposed by Andrew Dalke '''
    seen = set()
    for x in coll:
        if x not in seen:
            seen.add(x)
            yield x


def take(n, coll):
    new_coll = []
    for i in islice(coll, n):
        new_coll.append(i)
    return new_coll


def repeatedly(f, *args):
    result = f(*args)
    while True:
        yield result
        result = f(*args)


