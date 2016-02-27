from django.conf import settings
import importlib
from build_model import prepare_build_display_model

response_headers = ["boolean_crystallisation_outcome"]
descriptor_headers = ["reaction_temperature"]

visitorModules = {library:importlib.import_module(settings.STATS_MODEL_LIBS_DIR + "."+ library) for library in settings.STATS_MODEL_LIBS}


for modelVisitorLibrary, module in visitorModules.items():
    for modelVisitorTool in module.tools:
        print modelVisitorTool
        prepare_build_display_model(descriptor_headers, response_headers, modelVisitorLibrary, modelVisitorTool, "SingleSplitter", verbose=True)
