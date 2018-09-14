
import tellurium as te

# From Biotapestry: Export > Export to SBML
def convert_biotapestry_to_antimony(SBML_filename):
    return te.sbmlToAntimony(SBML_filename)

