import itertools

from chemlib.chemistry import Compound, Reaction
from chemlib.utils import reduce_list
import numpy

# TODO: TEMP IMPLEMENTATIONS, NEED TO CONSIDER ADDITIONAL METHODS/CALCULATIONS  
# INTERFACE?

class Combination(Reaction):
    '''
    Takes in the chemical symbol of the compounds as input,
    Then converts them into "Element" class objects which has all the properties,
    then tries to predict what will the products if they undergo
    a combination or synthesis reaction based on the element properties.
    '''

    def __init__(self, reactants=[Compound, Compound]):
        for n, compound in enumerate(reactants):
            if type(compound) is str:
                reactants[n] = Compound(compound)
        super().__init__(reactants, products=[])
        #self.balance()

    def is_complete(self):
        if not (self.products and self.products) : return False
        return True

    def _check_elemental_reactants(self, reactants_as_elements):
        for e in reactants_as_elements:
            if e.Group == 'Noble Gas': print("Noble gas present in list of elements: " + e.Element); return False
            if e.Valency == '|': print("Valancy not availible for: " + e.Element); return False
            else: return True

    def _get_oxidation_states(self, reactants_as_elements):
        oxidation_states = dict()
        for e in reactants_as_elements:
            oxidation_states[e.Element] = list(map(int, e.Oxidation.split('|')))
        return oxidation_states
    
    def _get_valid_oxidation_combinations(self, reactants_as_elements):
        '''
        When hydrogen and oxygen react with each other, they can form two main products:
        water (H2O) and hydrogen peroxide (H2O2).
        These reactions occur according to the following equations:

        2H2 + O2 -> 2H2O
        2H2 + 2O2 -> 2H2O2

        The first equation represents the formation of water from hydrogen and oxygen.
        In this reaction, the hydrogen atoms have an oxidation state of +1 and the oxygen atom
        has an oxidation state of -2. This reaction is exothermic, meaning that it releases heat.

        The second equation represents the formation of hydrogen peroxide from hydrogen and oxygen.
        In this reaction, the hydrogen atoms have an oxidation state of 0 and the oxygen atoms have an oxidation state of -1.
        This reaction is also exothermic.

        The electronegativities of hydrogen and oxygen play a role in determining the products of the reaction.
        Hydrogen has a lower electronegativity than oxygen,
        so it tends to have a positive oxidation state in compounds with oxygen.
        Oxygen, on the other hand, has a higher electronegativity, so it tends to have a negative oxidation state in
        compounds with hydrogen. This helps to explain why the oxidation states of hydrogen and oxygen in
        the equations above are consistent with their electronegativities.
        '''
        # get name of element with higher electronegetivity
        if reactants_as_elements[0].Electronegativity > reactants_as_elements[1].Electronegativity:
            element_with_higher_electronegetivity = reactants_as_elements[0].Element
        else:
            element_with_higher_electronegetivity = reactants_as_elements[1].Element

        # remove positive oxidation states for more electronegetive element
        # remove negative oxidation states for lesser electronegetive element
        # also remove 0 oxidation state
        oxidation_states_for_reactants = self._get_oxidation_states(reactants_as_elements)
        for element, oxidation_states in oxidation_states_for_reactants.items():
            for s in oxidation_states:
                if element_with_higher_electronegetivity == element and s > 0: oxidation_states.remove(s)
                if element_with_higher_electronegetivity != element and s < 0: oxidation_states.remove(s)

        oxidation_combinations = list([i for i in itertools.product(*oxidation_states_for_reactants.values())])
        return oxidation_combinations


    def complete_elemental_reaction(self):
        '''
        How to predict the products of a synthesis reaction between
        two elements looking at the properties of the reactants' atomic number,
        spatial configuration, valency,  electronegativity and all of the
        oxidation states?
        '''
        if not self.reactants or len(self.reactants) < 2: print("Need atleast 2 reactants to proceed!"); return self
        if self.is_complete(): print("Reaction Already complete"); return self
        reactants_as_elements = [r.elements[0] for r in self.reactants]
        reaction_eligible_for_reactants = self._check_elemental_reactants(reactants_as_elements)
        if not reaction_eligible_for_reactants: return self
        #element_valancy = list(map(int, str(e.Valency).split('|')))
        valid_oxidation_combinations = self._get_valid_oxidation_combinations(reactants_as_elements)

        
        # find all combinations from the array above

        print("Oxidation combination for {e1} x {e2}: {oxidation_matrix}".format(e1=reactants_as_elements[0].Element,
                                                                                 e2=reactants_as_elements[1].Element,
                                                                                 oxidation_matrix=valid_oxidation_combinations)
            )

        return self

if __name__ == "__main__":
    empty_r = Combination(reactants=[])
    complete_empty_r = empty_r.complete_elemental_reaction()
    print('completed_formula: ' + complete_empty_r.formula + '\n')

    r_with_noble_gas = Combination(reactants=['He'])
    complete_r_with_noble_gas = r_with_noble_gas.complete_elemental_reaction()
    print('completed_formula: ' + complete_r_with_noble_gas.formula + '\n')

    r_with_no_valancy = Combination(reactants=['Hg', 'F'])
    complete_r_with_no_valancy = r_with_no_valancy.complete_elemental_reaction()
    print('completed_formula: ' + complete_r_with_no_valancy.formula + '\n')

    glucose = Compound("C₆H₁₂O₆")
    print('elements: ' + str(set([e.Element for e in glucose.elements])) + '\n')

    methane = Compound('C2H5')
    from chemlib import Combustion
    combustion_r = Combustion(methane)
    print('completed_formula: ' + combustion_r.formula + '\n')

    r  = Combination(reactants=['H', 'O'])
    complete_r = r.complete_elemental_reaction()
    print('completed_formula: ' + r.formula + '\n')
