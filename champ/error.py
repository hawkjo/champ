def fail(message, return_code=1):
    # Instead of dumping stack traces back the user we provide them a readable message and a return code
    print(message)
    exit(return_code)
