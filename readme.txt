The source code can be found in the file "01Knapsack.cpp"

************************************************************************************************

An example of a toy problem instance can be found in the file "toyInstance.txt"

The first line of this file contains two integers: n and c, which denote the number of items and knapsack capacity respectively.
The next n lines contain two space separated integers, which denote the value and weight of an item.

************************************************************************************************

The file "01KnapsackParameters.txt" contains the parameters that the algorithm will use.
The first line contains the number t, which denotes the amount of seconds for which the algorithm will run.
The second line contains the integer w, which denotes the beam width.


************************************************************************************************

To execute a problem instance, first compile the source code (we assume a linux environment):

g++ -g -std=c++11 -O2 01Knapsack.cpp -o 01KnapsackExecutable

Now, run the executable on the toy problem instance and write the output to a file called "MCTSOutput.txt"

./01KnapsackExecutable < toyInstance.txt > MCTSOutput.txt

