#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <chrono>

using namespace __gnu_pbds;
using namespace std;

long long inf=1e18;
double eps=1e-9;

int width;

const int m=2; // for 01-knapsack, every variable has 2 choices: 0 (don't take) or 1 (take)

int n;
long long capacity;
vector<long long> weights;
vector<long long> values;
map<int, int> oldIdxToNewIdx;
int fractionPtr;
long long numerator;
long long denominator;

long long bestValue;
long long bestWeight;
vector<int> bestSolution;

long long currentWeight;
long long currentValue;
vector<int> currentSolution;

long long greedyWeight;
long long greedyValue;
vector<int> greedySolution;

double allowedTimeInSeconds;

random_device rd;

int randIntegralBetween(int lo, int hi)
{
    mt19937 gen(rd());
    uniform_int_distribution<> dis(lo,hi);
    return dis(gen);
}

double randDoubleZeroOne()
{
    return 1.0*randIntegralBetween(0, INT_MAX)/INT_MAX;
}

double randDouble(double lo, double hi)
{
    return lo+randDoubleZeroOne()*(hi-lo);
}


struct node
{
    vector<node*> children;
    double value;
    long long amountVisited;
    bool blocked;
    int numBlockedChildren;
    node* parent;
    node()
    {
        children.assign(m,NULL);
        value=inf;
        amountVisited=0;
        blocked=false;
        numBlockedChildren=0;
        parent=NULL;
    }
};
node* root;

void readParameters()
{
    ifstream inFile;
    inFile.open("01KnapsackParameters.txt");
    inFile >> allowedTimeInSeconds;
    inFile >> width;
    inFile.close();
}

void readProblemInstance() // from standard input
{
    scanf("%d %lld", &n, &capacity);
    vector< pair<double, pair<long long, pair<long long, int> > > > sorted_by_ratio; // <ratio, <val, <wei, i> > >
    int initial_n=n;
    for(int i=0; i<initial_n; i++)
    {
        long long wei, val;
        //scanf("%lld %lld", &wei, &val);
        scanf("%lld %lld", &val, &wei);
        if(wei>capacity)
        {
            n--;
            continue;
        }
        sorted_by_ratio.push_back(make_pair(1.0*val/wei,make_pair(val,make_pair(wei,i))));
    }
    sort(sorted_by_ratio.begin(),sorted_by_ratio.end());
    reverse(sorted_by_ratio.begin(),sorted_by_ratio.end()); // biggest val/wei earlier
    values.assign(n,-1);
    weights.assign(n,-1);
    //cerr << "n: " << n << endl;
    //cerr << "sortedByRatio.size() " << sorted_by_ratio.size() << endl;
    for(int i=0; i<n; i++)
    {
        values[i]=sorted_by_ratio[i].second.first;
        weights[i]=sorted_by_ratio[i].second.second.first;
        oldIdxToNewIdx[sorted_by_ratio[i].second.second.second]=i;
    }   
}

// update both current and greedy solution!
void makeChoice(int index, int i)
{
    currentSolution[index]=i;
    currentValue += i*values[index];
    currentWeight += i*weights[index];

    if(i==0)
    {
        /*if(0<=fractionPtr && fractionPtr<n)
        {
            greedySolution[fractionPtr]=0;
            greedyValue -= (numerator*values[fractionPtr])/denominator;
            greedyWeight -= (numerator*weights[fractionPtr])/denominator;
        }*/
        
        if(index < fractionPtr)
        {
            greedySolution[index]=0;
            greedyValue -= values[index];
            greedyWeight -= weights[index];
        }
    }
}

long long solutionCompleted() // function is only called if all elements are chosen (index >= n) or all remaining elements can be included (fractionPtr >= n)
{
    if(0<=fractionPtr && fractionPtr<n)
    {
        greedyValue -= (numerator*values[fractionPtr])/denominator;
    }
    /*cerr << "greedy value: " << greedyValue << endl;
    cerr << "and its solution: " << endl;
    for(int i=0; i<n; i++)
    {
        cerr << greedySolution[i] << " ";
    }
    cerr << endl;
    cerr << endl;*/

    if(greedyValue > bestValue)
    {
        bestValue=greedyValue;
        cerr << "new best value of: " << bestValue << endl;
        bestSolution=greedySolution;
        /*for(int i=0; i<bestSolution.size(); i++)
        {
            cerr << bestSolution[i] << " ";
        }
        cerr << endl;*/
    }
    return greedyValue;
}

int debCtr=0;
long long calcUpperBound(int index) // everything up until index (inclusive) has already been decided
{
    if(0<=fractionPtr && fractionPtr<n) // remove fractional item
    {
        greedyWeight-=(numerator*weights[fractionPtr])/denominator;
        greedyValue-=(numerator*values[fractionPtr])/denominator;
    }
    if(fractionPtr<=index)
    {
        fractionPtr=index+1;
    }
    while(fractionPtr<n && greedyWeight+weights[fractionPtr]<=capacity)
    {
        greedyWeight+=weights[fractionPtr];
        greedyValue+=values[fractionPtr];
        /*if(debCtr<5)
        {
            cerr << "values[" << fractionPtr << "]= " << values[fractionPtr] << endl;
            debCtr++;
        }*/
        greedySolution[fractionPtr]=1;
        fractionPtr++;
    }
    if(fractionPtr<n)
    {
        numerator=capacity-greedyWeight;
        denominator=weights[fractionPtr];
        greedyWeight=capacity;
        greedyValue+=(numerator*values[fractionPtr])/denominator;
    }
    return greedyValue;
}

long long completePartialSolution(int index)
{
    if(0<=fractionPtr && fractionPtr<n)
    {
        greedyValue -= (numerator*values[fractionPtr])/denominator;
        greedyWeight -= numerator;
        numerator=0;
        denominator=1;
    }
    for(int i=index; i<n; i++)
    {  
        if(greedyWeight+weights[i]<=capacity && greedySolution[i]==0)
        {
            greedyWeight+=weights[i];
            greedyValue+=values[i];
            greedySolution[i]=1;
        }
    }
    return solutionCompleted();
}

void cascadingBlock(node* node)
{
//    cerr << "root has " << root->numBlockedChildren << " blocked children" << endl;
//    if(node==root) cerr << "root will be blocked" << endl;
    if(node->blocked) return;
    node->blocked=true;
    if(node->parent != NULL)
    {
        node->parent->numBlockedChildren += 1;
        if(node->parent->numBlockedChildren==m) cascadingBlock(node->parent);
    }
}


int lastLevelCut=0;
vector<node*> lastLevel;
// -inf means the upper bound is too low
// -1 means some child returned -inf
long long buildSolution(int index, node* currentNode, bool firstTime=false)
{
    if(firstTime) return completePartialSolution(0);
    if(index>=n || fractionPtr>=n)
    {
        //cerr << "hierzo!" << endl;
        return solutionCompleted();
    }
    long long atMost = calcUpperBound(index-1);  
    if(atMost<=bestValue)
    {
        currentNode->value=inf;
        /*cerr << "index: " << index << endl;
        for(int i=0; i<=index; i++)
        {
            cerr << "sol[" << i << "]: " << greedySolution[i] << endl;
        }
        cerr << atMost << " " << bestValue << endl;*/
        return -inf;
    }

    if(index>=lastLevelCut)
    {
        vector< pair<double, int> > sortedValues;
        int amountEligible=0;
        for(int i=0; i<2; i++)
        {
            if(((i==0) || (i==1 && currentWeight+weights[index]<=capacity)) && (currentNode->children[i]==NULL || !currentNode->children[i]->blocked))
            {
                amountEligible++;
            } 
            else continue;
            if(currentNode->children[i] != NULL)
            {
                double value=currentNode->children[i]->value;
                sortedValues.push_back(make_pair(value,i));
            }
        }
        sort(sortedValues.begin(),sortedValues.end());
        reverse(sortedValues.begin(),sortedValues.end());

        vector<double> inverseRank(m,-1);
        for(int i=0; i<sortedValues.size(); i++)
        {
            inverseRank[sortedValues[i].second]=sortedValues.size()-i;
        }
        vector<double> childrenScores;
        double normalizationFactor=0;
        for(int i=0; i<2; i++)
        {
            if(!(((i==0) || (i==1 && currentWeight+weights[index]<=capacity)) && (currentNode->children[i]==NULL || !currentNode->children[i]->blocked)))
            {
                continue;
            } 
            if(currentNode->children[i] != NULL)
            {
                double score=inverseRank[i];
                childrenScores.push_back(score);
                normalizationFactor+=score;
            }
        }
        
        int scorePtr=0;
        double sumOfScores=0;
        pair<double, int> toMaximize=make_pair(-inf,-1);
        for(int i=0; i<2; i++)
        {
            if(!(((i==0) || (i==1 && currentWeight+weights[index]<=capacity)) && (currentNode->children[i]==NULL || !currentNode->children[i]->blocked)))
            {
                continue;
            } 
            if(currentNode->children[i] != NULL)
            {
                childrenScores[scorePtr]/=normalizationFactor;
                childrenScores[scorePtr]+=sqrt(2*log(currentNode->amountVisited)/currentNode->children[i]->amountVisited);
                toMaximize=max(toMaximize,make_pair(childrenScores[scorePtr],i));
                sumOfScores+=childrenScores[scorePtr];
                scorePtr++;
            }
        }
        scorePtr=0;
        vector<double> percentages(m);
        for(int i=0; i<m; i++)
        {
            if(!(((i==0) || (i==1 && currentWeight+weights[index]<=capacity)) && (currentNode->children[i]==NULL || !currentNode->children[i]->blocked)))
            {
                percentages[i]=0.0;
            } 
            else
            {
                if(currentNode->children[i] != NULL)
                {
                    percentages[i]=(1.0*childrenScores.size()/amountEligible)*childrenScores[scorePtr]/sumOfScores;
                    scorePtr++;
                }
                else
                {
                    percentages[i]=1.0/amountEligible;
                }
            }
            if(i-1>=0) percentages[i]+=percentages[i-1];
        }
        if(abs(1-percentages[m-1]) > eps)
        {
            //cerr << "should this happen?" << endl;
            return -inf; // this means no child can be explored    
        }
        percentages[m-1]=1.0; // avoid potential numerical errors
                

        double randNumber=randDoubleZeroOne();
        for(int i=0; i<m; i++)
        {
            if(randNumber<=percentages[i] && ((i==0 && randNumber != 0) || percentages[i-1]<randNumber))
            {
                if(currentNode->children[i]!=NULL) // choose child deterministically
                {
                    i=toMaximize.second;
                }
                makeChoice(index, i);   
                long long ret;
                if(currentNode->children[i] != NULL)
                {
                    ret=buildSolution(index+1, currentNode->children[i]);
                }
                else
                {
                    /*cerr << "entered!" << endl;
                    cerr << "chosen for item " << index << " is " << i << endl;*/
                    ret=completePartialSolution(index+1);
                    currentNode->children[i] = new node();
                    currentNode->children[i]->value=ret;
                    currentNode->children[i]->amountVisited=1;
                    currentNode->children[i]->parent=currentNode;
                }
                if(ret==-1)
                {
                    return ret;
                }
                if(ret==-inf)
                {
                    //cerr << "called from index " << index << endl;
                    cascadingBlock(currentNode->children[i]);
                    return -1;
                }
                currentNode->value += 1.0*ret/currentNode->amountVisited;
                currentNode->value *= (1.0*currentNode->amountVisited)/(currentNode->amountVisited+1);
                currentNode->amountVisited++;
                return ret;
            }
        }
    }
    else
    {
        vector< pair<double, int> > sortedValues;
        for(int i=0; i<2; i++)
        {
            if(currentNode->children[i] != NULL && !currentNode->children[i]->blocked)
            {
                double value=currentNode->children[i]->value;
                sortedValues.push_back(make_pair(value,i));
            }
        }
        sort(sortedValues.begin(),sortedValues.end());
        reverse(sortedValues.begin(),sortedValues.end());

        vector<double> inverseRank(m,-1);
        for(int i=0; i<sortedValues.size(); i++)
        {
            inverseRank[sortedValues[i].second]=sortedValues.size()-i;
        }
        vector<double> childrenScores;
        double normalizationFactor=0;
        for(int i=0; i<m; i++)
        {
            if(currentNode->children[i] != NULL  && !currentNode->children[i]->blocked)
            {
                double score=inverseRank[i];
                childrenScores.push_back(score);
                normalizationFactor+=score;
            }
        }
        
        pair<double, int> toMaximize=make_pair(-inf,-1);
        double sumOfScores=0;
        int scorePtr=0;
        for(int i=0; i<2; i++)
        {
            if(currentNode->children[i] != NULL  && !currentNode->children[i]->blocked)
            {
                childrenScores[scorePtr]/=normalizationFactor;
                childrenScores[scorePtr]+=sqrt(2*log(currentNode->amountVisited)/currentNode->children[i]->amountVisited);
                toMaximize=max(toMaximize,make_pair(childrenScores[scorePtr],i));
                sumOfScores+=childrenScores[scorePtr];
                scorePtr++;
            }
        }

        scorePtr=0;
        vector<double> percentages(m);
        for(int i=0; i<m; i++)
        {
            if(currentNode->children[i] != NULL && !currentNode->children[i]->blocked && currentWeight+weights[index]*i<=capacity)
            {
                percentages[i]=childrenScores[scorePtr]/sumOfScores;
                scorePtr++;
            }
            else
            {
                percentages[i]=0;
            }
            if(i-1>=0) percentages[i]+=percentages[i-1];
        }
        percentages[m-1]=1.0; // avoid potential numerical errors
        if(abs(1-percentages[m-1]) > eps)
        {
            //cerr << "should this happen?" << endl;
            return -inf; // this means no child can be explored
        } 

        double randNumber=randDoubleZeroOne();
        for(int i=0; i<m; i++)
        {
            if(randNumber<=percentages[i] && ((i==0 && randNumber !=0) || percentages[i-1]<randNumber))
            {
                i=toMaximize.second;
                makeChoice(index, i);
                long long ret=buildSolution(index+1, currentNode->children[i]);
                if(ret==-1)
                {
                    return ret;
                }
                if(ret==-inf)
                {
                    //cerr << "called from index: " << index << endl;
                    cascadingBlock(currentNode->children[i]);
                    return -1;
                }
                currentNode->value += 1.0*ret/currentNode->amountVisited;
                currentNode->value *= (1.0*currentNode->amountVisited)/(currentNode->amountVisited+1);
                currentNode->amountVisited++;
                return ret;
            }
        }
    }
}

long long generateSingleSample(bool firstTime=false)
{
    //cerr << "sample entered" << endl;
    fractionPtr=-2; // could be improved by not redoing this work

    greedyWeight=0;
    greedyValue=0;
    greedySolution.assign(n,0);
    
    currentWeight=0;
    currentValue=0;
    currentSolution.assign(n,0);
    if(root->blocked) return inf;
    return buildSolution(0, root, firstTime);
}

void performCut()
{
    priority_queue< pair<double, node*> > pq;
    for(node* parentPtr : lastLevel)
    {
        parentPtr->numBlockedChildren=m;
        for(int i=0; i<m; i++)
        {
            node* child=parentPtr->children[i];
            if(child==NULL)
            {
                continue;
            }
            child->blocked=true;
            pq.push(make_pair(child->value, child));
            if(pq.size()>width) pq.pop();
        }
    }
    vector<node*> newLastLevel;
    while(!pq.empty())
    {
        pair<double, node*> top=pq.top();
        pq.pop();
        node* nd=top.second;
        newLastLevel.push_back(nd);
        nd->blocked=false;
        nd->parent->numBlockedChildren -= 1;
    }
    for(node* parentPtr : lastLevel)
    {
        if(parentPtr->numBlockedChildren==m)
        {
//            cerr << "this was entered" << endl;
             cascadingBlock(parentPtr);
        }
    }
    lastLevel=newLastLevel;
//    cerr << "after perform cut, root has " << root->numBlockedChildren << " blocked children" << endl;
}

void searchLoop()
{
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    root = new node();
    lastLevel.push_back(root);
    bestWeight=0;
    bestValue=0;
    bestSolution.assign(n,0);
    generateSingleSample(true);
    while(true)
    {
//        cerr << "lastLevelCut: " << lastLevelCut << endl;
        generateSingleSample();
        if(root->blocked)
        {
//            cerr << "the rootwas blocked" << endl;
            break;
        }
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        double elapsedTime=chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0;
        if(elapsedTime>=allowedTimeInSeconds)
        {
//            cerr << "time is up!" << endl;
            return;
        }
//        cerr << "elapsedTime: " << elapsedTime << " allowedTimeInSeconds " << allowedTimeInSeconds << endl;
        if(elapsedTime>=(1.0*(lastLevelCut+1)/n)*allowedTimeInSeconds)
        {
            performCut();
            lastLevelCut++;
        }
    }
}

void precompute()
{
}

void debugVerifyIndeedCorrect()
{
    long long val=0;
    long long wei=0;
    for(int i=0; i<n; i++)
    {
        val+=values[i]*bestSolution[i];
        wei+=weights[i]*bestSolution[i];
    }
    if(wei>capacity || val != bestValue)
    {
        freopen("invalidCertificate.txt","w",stdout);
        cerr << "your certificate is invalid!" << endl;
    }
}

int main()
{
    bestValue=-inf;
    readParameters();
    readProblemInstance();
    precompute();
    /*for(int i=0; i<n; i++)
    {
        cerr << "item " << i << " has weight " << weights[i] << " and value " << values[i] << endl;
    }*/
    searchLoop();
    cout << "Best found solution has a value of: " << bestValue << endl;
    cerr << "best val: " << bestValue << endl;
    debugVerifyIndeedCorrect();
    return 0;
}
