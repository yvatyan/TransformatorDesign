#include <iostream>
#include <string>
using namespace std;

#define SQR(x) ((x)*(x))

enum material {
    Cu,
    Fe
}

struct Box {
    double height;
    double length;
    double width;
};

double metal_density (material metal) { // -> kg/m^3
    switch(metal) {
        case Cu :   return 8930.;
        case Fe :   return 7850.;
        default :   return -1;
    };
}

double metal_resistivity (material metal) { // -> Ohm*m [at 20C]
    switch(metal) {
        case Cu :   return 1.68e-8;
        case Fe :   return 9.71e-8;
        default :   return -1;
    };
}

double material_heat_capacity (material mat) { // -> J/(kg*K) 
     switch(mat) {
        case Cu :   return 385;
        case Fe :   return 444;
        default :   return -1;
    };
}

double wire_mass (double length, double diameter, material metal) { // m, m -> kg
    double V = SQR(diameter/2.)*M_PI*length;
    return V*metal_density(metal);
}

double heat (double mass, double temperature_diff, material mat) { // kg -> J
    return material_heat_capacity(mat)*temperature_diff*mass; 
}

double wire_resistance (double length, double diameter, material metal) { // m, m, -> Ohm
    double S = SQR(diameter/2.)*M_PI;
    return metal_resistivity(metal)*length/S; 
}

double coil_inductivity (double wireLength, double wireDiameter, Box coilSize, material coreMetal) { // m, m, m x m x m -> H

}

double coil_resistance (double wireLength, double wireDiameter, Box coilSize, double frequency, material coreMetal, material wireMetal) { // m, m, m x m x m, Hz -> Ohm
    return wire_resistance(wireLength, wireDiameter, wireMetal) + frequency*coil_inductivity(wireLength, wireDiameter, coilSize, coreMetal);
}

double coil_heat (double resistance, double current) { // Ohm, A -> J/sec
    return resistance*SQR(current);
}

int main() {

    

}
