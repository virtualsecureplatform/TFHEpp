# Tutorial

This directory contains the tutorial codes which demostrate how to use TFHEpp to resolve 2-party [Yao's Millionaires' problem](https://en.wikipedia.org/wiki/Yao%27s_Millionaires%27_problem). This is one of the simplest and most famous Secure Function Evaluation problem.

# How to build
When building TFHEpp, give `-DENABLE_TUTORIAL=ON` option to cmake.

# How to run
1. Run "client" and type the integer which represents client's wealth. "client" will produce encrypted input, the secret key and the gate key.
2. Run "cloud" and type the integer which represents cloud's wealth. "cloud" will do subtraction of two integer to compare integers and return the encrypted result.
3. Run "verify". This will decrypt the result and print it.
