#include <iostream>
#include <string>
using namespace std;

#define SQR(x) ((x)*(x))

namespace constant {
    const double pi = M_PI;
    const double u0 = 4*pi*1e-7;
}

enum material {
    Cu,
    Fe
}
enum packingType {
    quadratic,
    compact
}

struct Box {
    double height;
    double length;
    double width;
    double perimeter() const {
        return 4*(height + length + width);
    }
    double widthSizePerimeter() const {
        return 2*(height + width);
    }
    double lengthSizePerimeter() const {
        return 2*(height + length);
    }
    double basePerimeter() const {
        return 2*(length + width);
    }
    double area() const {
        2*(height*length + length*width + width*height);
    }
    double widthSideArea() const {
        return height*width;
    }
    double lengthSideArea() const {
        return height*length;
    }
    double baseArea() const {
        return length*width;
    }
    double volume() const {
        return height*length*width;
    }
};

struct CoilParams {
    double turns;
    double windingThickness;
}

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

double material_heating_capacity (material mat) { // -> J/(kg*K) 
     switch(mat) {
        case Cu :   return 385;
        case Fe :   return 444;
        default :   return -1;
    };
}

double material_permeability(material mat) {
    switch(mat) {
        case Cu :   return 1.26e-6;
        case Fe :   return 6.30e-3;
        default :   return -1;
    };
}

double wire_mass (double length, double diameter, material metal) { // m, m -> kg
    double V = SQR(diameter/2.)*constant::pi*length;
    return V*metal_density(metal);
}

double heating (double mass, double temperature1, double temperature2, material mat) { // kg -> J
    return material_heating_capacity(mat)*(temperature2 - temperature1)*mass; 
}

double wire_resistance (double length, double diameter, material metal) { // m, m, -> Ohm
    double S = SQR(diameter/2.)*constant_pi;
    return metal_resistivity(metal)*length/S; 
}

double winding (double wireLength, double wireDiameter, Box coreSize, packingType packing = quadratic) { // m, m, (m, m, m) -> count
    double k = coreSize.height/wireDiameter; 
}

double coil_inductivity (double wireLength, double wireDiameter, Box coreSize, material coreMetal) { // m, m, (m, m, m) -> H
    CoilParams params = winding(wireLength, wireDiameter, coreSize);
    double N = params.turns;
    double W = params.windingThickness;
    double R = W + sqrt(coreSize.baseArea() / constant::pi);
    return (material_permeability(coreMetal)*SQR(N)*SQR(R))/(2*constant::pi*(6*R + 9*coreSize.height + 10*W)); 
}

double coil_resistance (double wireLength, double wireDiameter, Box coilSize, double frequency, material coreMetal, material wireMetal) { // m, m, (m, m, m), Hz -> Ohm
    return wire_resistance(wireLength, wireDiameter, wireMetal) + frequency*coil_inductivity(wireLength, wireDiameter, coilSize, coreMetal);
}

double wire_heating (double resistance, double current) { // Ohm, A -> J/sec
    return resistance*SQR(current);
}

int main() {

    

}
