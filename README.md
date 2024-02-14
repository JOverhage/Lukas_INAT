# Lukas_INAT
Lukas '78 asset pricing model with inattention

# 1st Idea:
Unclear if this sort of model conforms to aggregate linearity, approximate aggregation etc.
For now, krussel-smith ish idea: Guess that price is just a fct of aggregate state TODAY 
Under this, can guess two prices, solve VFIs for this, then simulate & check if under these prices, asset mkt always clears
If not, update guess for price and repeat
