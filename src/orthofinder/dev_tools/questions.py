import inspect 


def question(message: str, 
                  counter_ctrl: list[str] = [""],
                  question_no: list[int] = [0]) -> None:
    
    f = inspect.currentframe()
    info = inspect.getframeinfo(f.f_back)
    question_message = f"{info.filename}:{info.function}:{info.lineno+1}: - {message}"
    
    if question_message not in counter_ctrl:
        question_no[0] += 1
        counter_ctrl[0] = question_message
        print()
        print('*'*50)
        print(f">>> {question_no[0]}. QUESTION:")
        print(question_message)
        print('*'*50)


