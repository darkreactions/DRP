def get_model(model_type):

  if model_type=="random forest":
    from sklearn.ensemble import RandomForestClassifier as model
    descriptors = {"n_estimators":500, "criterion":"entropy", "n_jobs":-1}

  elif model_type=="linear regression": #TODO: "Cannot perform reduce with flexible type"
    from sklearn.linear_model import LinearRegression as model
    descriptors = dict()

  elif model_type=="svc":
    from sklearn.svm import SVC as model
    descriptors = {"C":1, "kernel":"linear"}

  elif model_type=="knn":
    from sklearn.neighbors import KNeighborsClassifier as model
    descriptors = {"p":3, "n_neighbors":1, "weights":"distance"}

  else:
    raise Exception("Model model_type '{}' unknown by get_model".format(model_type))

  return model(**descriptors), descriptors

