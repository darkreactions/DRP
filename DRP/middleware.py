"""Middleware classes for handling specific cases in DRP."""
from django.shortcuts import render
from chemspipy.errors import ChemSpiPyServerError
import logging

logger = logging.getLogger(__name__)

class ChemspiderErrorMiddleware(object):

    """
    Handling for Chemspider Server Errors.

    From time to time, the chemspider server will have a bug which is beyond 
    our control, in order to cope with this, this class will emit a more
    helpful 500 page than the standard.
    """

    def processException(request, exception):
        """The handler for the exception."""
        if isinstance(exception, ChemSpiPyServerError):
            logger.error(exception.message)
            render(request, "chemspider_500.html")
        else:
            return None
