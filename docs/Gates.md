# Gates in TFHEpp

On this page, some logical gates supported in TFHEpp are introduced. Most of gates are defined in [gate.hpp](../include/gate.hpp). 
As you may know, TFHE is known as one kind of FHE which is suitable for logical gate operations, so you may want to start with using these operations for your application. 

## 2-input-1-output gates

Becuase all of the gates classified into this category(AND,XOR,OR,NAND,XNOR,NOR,ANDYN,ANDNY,ORYN,ORNY) has same interface, only the interface for NAND gate is described here. 
As you can easily find, there are two definition of HomNAND. 

```
template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ, class iksP = lvl10param>
void HomNAND(TLWE<typename iksP::targetP> &res, const TLWE<typename brP::domainP> &ca, const TLWE<typename brP::domainP> &cb,
                    const EvalKey &ek)

template <class iksP = lvl10param, class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ>
void HomNAND(TLWE<typename brP::targetP> &res, const TLWE<typename iksP::domainP> &ca, const TLWE<typename iksP::domainP> &cb,
             const EvalKey &ek)
```

The difference between these two definition is that the order of Identity Key Switching and Blind Rotate. This will effect the noise growth of the operation. If you are not sure about the TFHE's theory, the lower one is recommended because it causes lower error growth. 
There are 4 inputs for this function, `res,ca,cb,ek`. 
`res` is  the ciphertext to store the result. `ca` and `cb` are input ciphertexts. `ek` is the evaluation key. 
You may wonder about the existence of the gates followed by `YN` or `NY` like `ANDYN`. These gate will negate one of thire inputs. For example, `ANDYN` is equivalent to `~a&b`. This kind of gates will be used in the synthesis tools like [Yosys](https://github.com/YosysHQ/yosys). 