# LICENSE_GOES_HERE
from multiprocessing.managers import BaseManager
import traceback
from multiprocessing.pool import Pool
import multiprocessing


# Shortcut to multiprocessing's logger
def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)


class MyManager(BaseManager):
    pass


class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable
        return

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise

        # It was fine, give a normal answer
        return result
    pass

class LoggingPool(Pool):
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)
    def map(self,func,iterable):
        return Pool.map(self, LogExceptions(func), iterable)
def go():
    print(1)
    raise Exception()
    print(2)


if __name__ == '__main__':
    multiprocessing.log_to_stderr()
    p = LoggingPool(processes=1)
    
    p.apply_async(go)
    p.close()
    p.join()
    pass