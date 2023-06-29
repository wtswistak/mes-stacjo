#include <iostream>
#include <cmath>

using namespace std;

struct GlobalData {
    double Tot;   // temperatura otoczenia
    double alfa;  // warunek brzegowy konwekcji
    double q;     // gestosc strumienia
    double k;     // wspolczynnik przewodzenia ciepla
    double L;     // dlugosc preta
    double S;
    int nN;       // l wezlow
    int nE;       // l elementow
};

struct Node {
    double x;
    int BC;       // 1 konwe 2 strumien ciepla
    double t;     // temp na wezle
};

struct Element {
    double k; // wspolczynnik przewodzenia ciepla
    int id[2];       // E1-N1,N2-id[1,2]
    double Le;       // dlugosc elemn
    double** Hlocal; // lokalna macierz H
    double* Plocal;
    Node* nodes;     // wskaznik do tablicy węzłów
};

struct SiatkaMes {
    Element* elements;
    Node* nodes;
};

struct SOE {
    double** Hglobal;
    double* Pglobal;
    double* T;
};

Node nodes[] = {
    {0, 1, 0},
    {2.5, 0, 0},
    {5, 2, 0}
};

void calcMacierzWektorLokalny(Element& element, GlobalData data) {
    element.Le = element.nodes[1].x - element.nodes[0].x;
    double C = data.S * data.k / element.Le;

    element.Hlocal = new double* [2];
    element.Plocal = new double[2];
    element.Plocal[0] = 0;
    element.Plocal[1] = 0;

    element.Hlocal[0] = new double[2];
    element.Hlocal[1] = new double[2];
    element.Hlocal[0][0] = C;
    element.Hlocal[0][1] = -C;
    element.Hlocal[1][0] = -C;
    element.Hlocal[1][1] = C;

    if (element.nodes[0].BC == 1) {
        element.Plocal[0] = data.q * data.S;
    }
    else if (element.nodes[0].BC == 2) {
        element.Plocal[0] = -data.alfa * data.S * data.Tot;
        element.Hlocal[0][0] += data.alfa * data.S;
    }

    if (element.nodes[1].BC == 1) {
        element.Plocal[1] = data.q * data.S;
    }
    else if (element.nodes[1].BC == 2) {
        element.Plocal[1] = -data.alfa * data.S * data.Tot;
        element.Hlocal[1][1] += data.alfa * data.S;
    }
}

void calcMacierzWektorGlobalny(SiatkaMes& siatka, GlobalData data, SOE& soe) {
    soe.Hglobal = new double* [data.nN];
    for (int i = 0; i < data.nN; i++) {
        soe.Hglobal[i] = new double[data.nN];
        for (int j = 0; j < data.nN; j++) {
            soe.Hglobal[i][j] = 0.0;
        }
    }

    for (int i = 0; i < data.nE; i++) {
        Element& element = siatka.elements[i];
        soe.Hglobal[element.id[0] - 1][element.id[0] - 1] += element.Hlocal[0][0];
        soe.Hglobal[element.id[0] - 1][element.id[1] - 1] += element.Hlocal[0][1];
        soe.Hglobal[element.id[1] - 1][element.id[0] - 1] += element.Hlocal[1][0];
        soe.Hglobal[element.id[1] - 1][element.id[1] - 1] += element.Hlocal[1][1];
    }

    soe.Pglobal = new double[data.nN];
    for (int i = 0; i < data.nN; i++) {
        soe.Pglobal[i] = 0;
    }

    for (int i = 0; i < data.nE; i++) {
        Element& element = siatka.elements[i];
        soe.Pglobal[element.id[0] - 1] += element.Plocal[0];
        soe.Pglobal[element.id[1] - 1] += element.Plocal[1];
    }
}

void gauss(SOE& soe, int n) {
    double** tab = soe.Hglobal;
    double* rozw = soe.Pglobal;
    double* temp = new double[n];

    for (int j = 0; j < n - 1; j++) {
        for (int i = j + 1; i < n; i++) {
            double mnoznik = tab[i][j] / tab[j][j];
            for (int k = j; k < n; k++) {
                tab[i][k] = tab[i][k] - tab[j][k] * mnoznik;
            }
            rozw[i] = rozw[i] - rozw[j] * mnoznik;
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        temp[i] = rozw[i];
        for (int j = i + 1; j < n; j++) {
            temp[i] -= tab[i][j] * temp[j];
        }
        temp[i] /= tab[i][i];
    }

    for (int i = 0; i < n; i++) {
        soe.T[i] = temp[i];
    }

    delete[] temp;
}

int main() {
    GlobalData data = {
        400,   // tempOtoczenia
        10,    // alfa
        -150,  // q
        50,    // k
        5,     // L
        2,     // S
        3,     // nN
        2      // nE
    };
    SiatkaMes siatka;
    SOE soe;

    siatka.nodes = new Node[data.nN];
    for (int i = 0; i < data.nN; i++) {
        siatka.nodes[i] = nodes[i];
    }

    siatka.elements = new Element[data.nE];
    for (int i = 0; i < data.nE; i++) {
        Element& element = siatka.elements[i];
        element.k = data.k;
        element.id[0] = i + 1;
        element.id[1] = i + 2;

        element.nodes = new Node[2];
        element.nodes[0] = siatka.nodes[i];
        element.nodes[1] = siatka.nodes[i + 1];
    }

    for (int i = 0; i < data.nE; i++) {
        Element& element = siatka.elements[i];
        calcMacierzWektorLokalny(element, data);
    }

    calcMacierzWektorGlobalny(siatka, data, soe);
    cout << "Macierz globalna H:" << endl;
    for (int i = 0; i < data.nN; i++) {
        for (int j = 0; j < data.nN; j++) {
            cout << soe.Hglobal[i][j] << " ";
        }
        cout << endl;
    }
    cout << "\nGlobalny wektor P:" << endl;
    for (int i = 0; i < data.nN; i++) {
        cout << soe.Pglobal[i] << " " << endl;
    }

    soe.T = new double[data.nN];
    gauss(soe, data.nN);

    cout << "\nTemperatury:" << endl;
    for (int i = 0; i < data.nN; i++) {
        cout << "Temperatura " << i + 1 << ": " << abs(soe.T[i]) << endl;
    }


    for (int i = 0; i < data.nN; i++) {
        delete[] soe.Hglobal[i];
    }
    delete[] soe.Hglobal;

    for (int i = 0; i < data.nE; i++) {
        Element& element = siatka.elements[i];
        delete[] element.Hlocal;
        delete[] element.Plocal;
        delete[] element.nodes;
    }
    delete[] siatka.elements;
    delete[] siatka.nodes;


    return 0;
}
