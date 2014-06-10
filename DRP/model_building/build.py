import model_methods
import time

if __name__ == "__main__":
	model_name = "{0}.model".format(time.time())
	model_methods.gen_model(model_name)
	
