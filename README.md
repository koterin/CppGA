# C++ Genetic Algorithms
AI-based experimental data approximation curve selection.

This project lets you set some experimental data (in a form of two-dimensional graph) and select an approximation curve in a symbolic form - for example, for something like that

![ExpData](images/expData.jpg)

The output formula would be v(t)=((((t+0.2)^(8.1)*t+t)^(4.28)*t+t)^2+1.01)^(-0.68) with the final accuracy of 90%

To process symbolic operation an open lib is used.

## What is genetic algorithms

[Wikipedia](https://en.wikipedia.org/wiki/Genetic_algorithm) will probably tell it better than me, but I'll give an overview:
We know that our genes have some repetitive data that forms into us, individuals, and may be our sister and brothers - closest to us in one population. Something very similar can be extrapolated to the algorithms - if we call one formula an individual, group of similar formulas - a population, and every variable, operator and operand inside of it - a gene. Let's take a look at the set of actions to process our algorithm:
1. Create first population - for example, y = x + 1
2. Populate it with similar expressions - y = 2 * x + 2, y = 0.5 * x + 3, y = x + 4, ...
3. Randomly choose MOM & DAD and crossover them - let's take part 2 * x from MOM and + 3 from DAD. Now we get a KID: y = 2 * x + 3
4. Rate the KID's goodness by some mark - for example, how well he desribes the input experimental data. The function to rate that is called fitness function
5. Put the KID back into the population and repeat from the step 3 instead of the worst individual in population
6. After each loop you will be getting better and better population, until the final KID meets you criteria

After getting familiar with this algorithm, you can make it spicier:
- Add some mutation: what if after MOM and DAD had some fun night a gene like ^2 suddenly appeared in a KID? It's a nice method to warm up a population that was dealing with too much incest and got all very similar
- Set some rules for crossover: what if you're dealing with really big formulas and you'd like to slice some random part from MOM and random part from DAD?
- Play with the survival rules: may be the KID should take its place in the population only if he is better than the worst individual?

## What is done here

This project includes the following logic:
1. The initial population is set from the one starting formula, which maybe really complex or not so much (this option helps to make iterative runs - you can put your output formula from the previous run into the initial in the next run)
2. Each individual goes into solo cycle of numerical GA: the symbolic structure stays the same, for example, y = 2 * x + 3, but the coefficients may ruin the whole thing, so the same way working algorithm selects some numbers to turn the individual into something like y = 2.5 * x + 1.4
3. The fitness function calculates the goodness of each individual (here I use the normalized standart deviation turned into percentage rate, where 100% is prefectly fitted individual)
4. The crossover happens: MOM and DAD can't be the same individual, what part to take from them and from where is decided randomly, no zero or primitive KID is allowed
5. Some mutation included: if the probability (which can be set manually) is less than N, then mutation happens
6. The termination criteria is checked
7. If it hasn't been met, the KID goes back into the pit

## How to run

Just set all the neccesary input data and run the Source.cpp file:

Go for one-variable formula selection in the **source** branch or multi-variable version in the **multi_symb** branch.
You can add new input data file with lines like this
```
datafiles.push_back("Data\\angles\\1.txt");
```
Variables to be considered:
```
	Variables.push_back(v);
	Variables.push_back(t);
	Variables.push_back(h);
```
And their values:
```
	hValues.push_back(1);
	hValues.push_back(0.5);
	hValues.push_back(0.25);
	hValues.push_back(0.1);
```
Mind that complicated calculation REALLY TAKE TIME. The one given in example was running for 2 hours on i510, the ones with multiple variables - for days
