__copyright__ = 'Copyright (C) 2018 ICTP'
__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__credits__ = ["Stefano Piani"]


class Repeater(object):
    """
    A Repeater is substantially a wrapepr around a function.
    Its aim is to avoid to compute a function several times if the function
    is called with the same argument
    """
    def __init__(self, function):
        self._function = function
        self._previous_value = None
        self._previous_input = None
        self._already_used = False

    def __call__(self, x):
        if not self._already_used or self._previous_input != x:
            self._previous_input = x
            self._previous_value = self._function(x)
            self._already_used = True

        return self._previous_value
