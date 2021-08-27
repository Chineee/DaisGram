#include <iostream>

using namespace std;

class sos {

    public:

        int a;
        int b;

        sos() {
            a = 666;
            b = 333;
        }

        sos& operator=(const sos& rhs) {
            // codice

            this->a = rhs.a;
            this->b = rhs.b;

            return *this;
        }
};

int main() {

   // sos* a = (sos*) malloc(10*sizeof(sos));
   sos a;
   sos b;
   sos c;
   sos d;

   (cout << a.a) << " " << a.b << endl;



   b.a = 1;
   b.b = 2;

   a = b = c = d;  
   a = (b = ( c = d));

   cout << a.a << " " << a.b << endl;



    return 0;
}