
#include "solProducer.hpp"



int main() {
    int N = 3000;
    CONSTS con1 = CONSTS();
    con1.setN(N);
    con1.initDir();
    con1.setDk();
    producer prd = producer(con1);
    double s = 1;
    prd.print123(s);




    return 0;
}
