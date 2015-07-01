class PonyToolsError(Exception):
    pass

class TriAllelicError(PonyToolsError):
    def __init__(self,expr,message=''):
        self.expr = expr
        self.message = message
