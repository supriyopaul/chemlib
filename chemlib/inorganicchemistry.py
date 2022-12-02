from chemlib.chemistry import Compound, Reaction
from chemlib.utils import reduce_list
import numpy as np

# TODO: TEMP IMPLEMENTATIONS, NEED TO CONSIDER ADDITIONAL METHODS/CALCULATIONS  
# INTERFACE?

class Combination(Reaction):

    def __init__(self, reactants=[Compound, Compound]):
        for n, compound in enumerate(reactants):
            if type(compound) is str:
                reactants[n] = Compound(compound)
        super().__init__(reactants, products=[])
        #self.balance()

    def is_complete(self):
        if not (self.products and self.products) : return False
        return True

    def complete_elemental_reaction(self):
        if self.is_complete(): print("Reaction Already complete"); return self
        reactants_as_elements = [r.elements[0] for r in self.reactants]
        for e in reactants_as_elements:
            element_valancy = e.Valency
            if element_valancy == '|': print("Valancy not availible for: " + e.Element); return self
            try:
                element_valancy = list(map(int, str(e.Valency).split('|')))
            except Exception as exp:
                import pdb; pdb.set_trace()
            if element_valancy == [0]: print("Noble gas present in list of elements: " + e.Element); return self
        
        oxidation_states = [list(map(int, e.Oxidation.split('|'))) for e in reactants_as_elements]
        print(oxidation_states)
        return self

if __name__ == "__main__":
    r_with_noble_gas = Combination(reactants=['He'])
    complete_r_with_noble_gas = r_with_noble_gas.complete_elemental_reaction()
    print('completed_formula: ' + complete_r_with_noble_gas.formula)

    r_with_no_valancy = Combination(reactants=['Hg', 'F'])
    complete_r_with_no_valancy = r_with_no_valancy.complete_elemental_reaction()
    print('completed_formula: ' + complete_r_with_no_valancy.formula)

    r  = Combination(reactants=['H', 'O'])
    complete_r = r.complete_elemental_reaction()
    print('completed_formula: ' + r.formula)
