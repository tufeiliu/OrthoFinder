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


# def question_func(message: str, 
#                   message_ctrl: list[tuple] = [],
#                   counter_ctrl: list[str] = [""],
#                   question_no: list[int] = [0]):

#     question_func.caller += 1
#     if (question_func.caller, message) not in message_ctrl:
#         message_ctrl.append((question_func.caller, message))
    
#     def question_decor(func):

#         @functools.wraps(func)
#         def wrapper(*args, **kwargs):
#             _, lineno = inspect.getsourcelines(func)
#             question_message = f"{inspect.getsourcefile(func)}:{func.__name__}:{lineno} - {message}"
#             if len(message_ctrl) != 0:
#                 for no, mes in message_ctrl:
#                     if (mes in question_message) and (question_message not in counter_ctrl):
#                         counter_ctrl.append(question_message)
#                         print()
#                         print('*'*50)
#                         print(f">>> {no}. FUNC QUESTION:")
#                         print(question_message)
#                         print('*'*50)

#             result = func(*args, **kwargs)
#             return result
#         return wrapper
#     return question_decor
# question_func.caller = 0
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
        
