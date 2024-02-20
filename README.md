# Lukas_INAT
Project to learn about how the Lukas '78 model behaves if inattention is added

# Starting point:
Unclear if this sort of model conforms to aggregate linearity, approximate aggregation etc.
For now, krussel-smith ish idea: Guess that price is just a fct of aggregate state TODAY 
Under this, can guess two prices, solve VFIs for this, then simulate & check if under these prices, asset mkt always clears
If not, update guess for price and repeat

# Check:
state-specific prices obtained in first step would be consistent if at the prices, full asset tree is demanded regardless of current history. On the other hand, if history matters there would be deviations in $A^s - A^d$ depending on how the economy has moved. Ex ante it seems reasonable these deviations would exist (more and more households update, e.g. more people are aware that state is bad if the economy has been in a bad state for 10 periods vs if it had been there only for one)

# Outcome:
Main Result was that not necessarily clear a bisection algorithm works, regular $p \uparrow \Rightarrow A^d \downarrow$ relation doesn't necessarily (monotonically?) persist! This might be resolved if the history is simply added as a state variable.
