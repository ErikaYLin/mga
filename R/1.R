# Set class for `mga` output
methods::setOldClass("mga")


new.mga <- methods::setClass("mga", contains = "list", representation = representation(info = "list"), prototype = methods::prototype(list(), info = list()))
