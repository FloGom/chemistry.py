class Element(object):
    def __init__(self, options):
        props = {
            "number"            : int,
            "symbol"            : str,
            "name"              : str,
            "mass"              : float,
            "exact_mass"        : float,
            "ionization"        : float,
            "electron_affinity" : float,
            "electronegativity" : float,
            "radius_vdw"        : float,
            "radius_covalent"   : float,
            "boiling_point"     : float,
            "melting_point"     : float,
            "block"             : str,
            "period"            : int,
            "group"             : str,
            "family"            : str,
            "electron_config"   : str
        }
        
        for prop in options:
            if prop not in props:
                raise AttributeError("Unknown property " + str(prop))
        
        for prop in props:
            if prop in options:
                setattr(self, prop, props[prop](options[prop]))
            else:
                setattr(self, prop, None)
                
        # more properties to come
    
    def __repr__(self):
        return self.name or self.symbol or ""

    __str__ = __repr__
    
    
