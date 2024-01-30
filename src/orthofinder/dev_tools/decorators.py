import functools
import time
import inspect

def timeit(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        time_elapsed = end - start
        print()
        print('***Function running time***')
        print(f"{func.__name__}: {time_elapsed:5f}s", end='\n'*2)
        return result
    return wrapper


def question_func(message):
    def question_decor(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            wrapper.calls += 1
            print()
            print('*'*50)
            print(f">>> {wrapper.calls}. FUNC QUESTION:")
            _, lineno = inspect.getsourcelines(func)
            print(f"{inspect.getsourcefile(func)}:{func.__name__}:{lineno} - {message}")
            print('*'*50)
            
            result = func(*args, **kwargs)
            return result
        wrapper.calls = 0
        return wrapper
    return question_decor

# def question_counter(func):
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         wrapper.calls += 1
#         if wrapper.calls < 2:
#             print()
#             print('*'*50)
#             print(f">>>CODE QUESTION {wrapper.calls}:")
#             result = func(*args, **kwargs)
#             print('*'*50, end='\n'*2)
#             return result

#     wrapper.calls = 0
#     return wrapper

# class Counter:
#     def __init__(self, func):
#         self.func = func
#         self.count = count

#     def __call__(self, *args, **kwargs):
#         self.count += 1
#         return self.func(*args, **kwargs)
        
