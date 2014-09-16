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
    
    
import xml.etree.cElementTree as ET

tree = ET.ElementTree(file="bodr-10/elements/elements.xml")
root = tree.getroot()

keymap = {
    "bo:atomicNumber"             : ("number",            "text"),
    "bo:symbol"                   : ("symbol",            "value"),
    "bo:name"                     : ("name",              "value"),
    "bo:mass"                     : ("mass",              "text"),
    "bo:exactMass"                : ("exact_mass",        "text"),
    "bo:radiusCovalent"           : ("radius_covalent",   "text"),
    "bo:radiusVDW"                : ("radius_vdw",        "text"),
    "bo:ionization"               : ("ionization",        "text"),
    "bo:electronAffinity"         : ("electron_affinity", "text"),
    "bo:electronegativityPauling" : ("electronegativity", "text"),
    "bo:boilingpoint"             : ("boiling_point",     "text"),
    "bo:meltingpoint"             : ("melting_point",     "text"),
    "bo:group"                    : ("group",             "text"),
    "bo:periodTableBlock"         : ("block",             "text"),
    "bo:period"                   : ("period",            "text"),
    "bo:family"                   : ("family",            "text"),
    "bo:electronicConfiguration"  : ("electron_config",   "text"),
}

for element in root[2:]:
    properties = {}
    for property in element:
        key = property.attrib["dictRef"]
        if key not in keymap:
            continue
        key = keymap[key]
        if key[1] == "text":
            properties[key[0]] = property.text
        else:
            properties[key[0]] = property.attrib[key[1]]
    
    globals()[properties["symbol"]] = Element(properties)

