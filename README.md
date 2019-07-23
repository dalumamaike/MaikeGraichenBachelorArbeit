# MaikeGraichenBachelorArbeit
Dieses Repository beinhaltet den Code für die Bachelorarbeit "Modellierung der Diffusion durch eine Membran in dreidimensionaler Auflösung" von Maike Graichen.

Bestandteil sind folgende Dateien:

Modell1Generator.vrlp
Modell1.ugx
Modell2.ugx
simulationDiffusion.lua

In den weiteren Abschnitten wird die Nutzung dieser Dateien beschrieben.

# VRL-Studio
Die Datei Modell1Generator.vrlp wurde mit VRL-Studio erstellt und muss auch mit diesem ausgeführt werden.
Weitere Informationen zu VRL-Studio sind hier zu finden:
https://vrl-studio.mihosoft.eu

# ProMesh
Das Programm wurde zur Erstellung der Modell1.ugx und Modell2.ugx Dateien genutzt.
Weitere Informationen sind hier zu finden:
http://www.promesh3d.com

# UG4
Weitere Informationen zu ug4 können sind hier zu finden:
https://github.com/UG4/ugcore

Based on UG4 (www.ug4.org/license)

The following bibliography is recommended for citation and must be preserved in all covered files: "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively parallel geometric multigrid solver on hierarchically distributed grids. Computing and visualization in science 16, 4 (2013), 151-164" "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel flexible software system for simulating pde based models on high performance computers. Computing and visualization in science 16, 4 (2013), 165-179"

# Ausführen des Skripts
Nach der erfolgreichen Kompilierung von ug4, kann das 'simulationDiffusion.lua'-Skript beispielsweise folgendermaßen ausgeführt werden:

cd ug4/runs 
ugshell -ex Examples/simulationDiffusion.lua

Mögliche weitere Parameter können folgende sein:

  -grid (Modell1 oder Modell2)
  -numRefs (Anzahl der Verfeinerungen)
  -dif (Werte auf der Hauptdiagonalen für den Diffusionstensor)
