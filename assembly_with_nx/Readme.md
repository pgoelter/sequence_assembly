# Bioinformatik - Übung 2
###Git Repository: https://github.com/pgoelter/sequence_assembly/tree/main/assembly_with_nx
# Aufgabe 1
### Files und Ordner
````
assembly_with_nx
|
|----data             // Enthält Beispieldateien mit Fragmenten zum Ausführen
|      |--frag.dat
|      |--frag1.dat
|      |--frag2.dat
|
|----graph            // Enthält alle Implementierungen bezüglich des overlap graphen
|      |--__init__.py // Implementierung für Knoten, Kanten und Graph
|      |--utils.py    // Implementierung von diversen Hilfsmethoden, wie z.B. String Überlappung
|      
|----assembler.py     // Einstiegspunkt. Implementierung des Parsing der Konsolenargumente.
|----Readme.md        // What you see right now!
|----environment.yml  // Alle dependencies in Form eines .yml files, mit welchem direkt ein Anaconda environment erstellt werden kann.
````
### Installation
Damit das Tool problemlos funktioniert muss die Software Graphviz installiert werden.  
Diese kann hier von der Website bezogen werden: https://graphviz.org/download/

1. Die Installation erfordert eine bestehende Anaconda Installation. Hierzu kann entweder Anaconda oder alternativ Miniconda verwendet werden um die virtuelle Umgebung zu erstellen.  
Anaconda Link: https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html  
Miniconda Link: https://docs.conda.io/en/latest/miniconda.html
   

2. Mithilfe des im Projektverzeichnis liegenden file **environment.yml** kann die Umgebung erstellt werden. Dazu in einer Konsole (z.B. Anaconda Prompt) in das Projektverzeichnis wechseln und dort den folgenden Befehl ausführen:
    ````bash
    conda env export > environment.yml
    ````  

3. Nach erfolgreicher Installation wurde die Umgebung mit Namen **bioinfo** installiert. Diese kann nun mit folgendem Befehl aktiviert werden:
    ````bash
    conda activate bioinfo
    ````
4. Installation abgeschlossen.

### Ausführen via IDE oder Konsole
Um das Tool auszuführen wird vorausgesetzt, dass die im Installationsschritt aufgeführte virtuelle Umgebung installiert wurde.  

#### Konsole
1. Umgebung aktivieren: 
    ````bash
    conda activate bioinfo
    ````
2. Ins Projektverzeichnis wechseln (../assembly_with_nx/)
3. Tool über Konsolenbefehl starten:
    ````bash
    python assembler.py "Pfad/zu/.dat/Datei"
    ````
4. Die Ausgabe erfolgt über die Konsole.
5. Für weitere Optionen kann die Hilfe aufgerufen werden:
    ````bash
    python assembler.py --help
    ````
   
#### IDE
1. Anaconda Umgebung in IDE einbinden (z.B. Pycharm)
2. Eintrittspunkt ist **assembly_with_nx/assembler.py**
3. Dort können entsprechend Pfad zur Datei, etc.. zum Ausprobieren auch von Hand eingetragen werden.

#### Starten mit der vorgefertigten .exe
Im Projektverzeichnis liegt ein eigenständiges executable (**assembler.exe**). Diese Datei kann analog wie im Schritt **Konsole** beschrieben über eine beliebige Konsole genutzt werden.  
**Hiermit kann das Tool auch ohne vorige Installation der virtuellen Umgebung ausgeführt werden.**
````bash
assembler.exe "data/frag1.dat"
````

**Mit weiterem Befehl können alle zwischenschritte in welchen der graph bearbeitet werden als pdf durch die Graphviz Software ausgegeben werden:**
````bash
assembler.exe --print_graph "data/frag1.dat"
````

**Oder auch nur der resultierende Graph. Hierbei handelt es sich bei erfolgreicher Zusammensetzung um einen einzelnen Knoten. Falls die Sequenz nicht zusammensetzbar ist sind die übrigen Knoten zu sehen:**
````bash
assembler.exe --print_only_result "data/frag1.dat"
````

** Alle mit der .exe Datei ausgeführten Befehle können analog auch mit python ausgeführt werden indem 'assembler.exe' durch 'python assembler.py' ersetzt wird.
#### Verfügbare Optionen
````bash
assembler.exe --help
usage: assembler.exe [-h] [--print_graphs] [--print_only_result] [--assemble_hamilton] [--assemble_greedy] path

Assemble DNA fragments from a *.frag file.

positional arguments:
  path                 The path containing the DNA fragments (reads).

optional arguments:
  -h, --help           show this help message and exit
  --print_graphs       Print all graphs in between each merge of two nodes.
  --print_only_result  Prints only the resulting graph. Should be a single node if everything worked.
  --assemble_hamilton  NOTE: CURRENTLY NOT WORKING! Todos: Calculate orientation; Updating the graph after
                       finding the hamilton path; Assembles the fragments by building the overlap graph, finding a
                       hamilton path with max summed up weight. Then merges all nodes of the path together.
  --assemble_greedy    Assembles the fragments by building the overlap graph and merging the nodes afterward by
                       picking the edges with the biggest weight and breaking ties arbitrarily.
````
# Aufgabe 2
**Anmerkung:** Nicht bearbeitet, da nur für 2er Gruppen vorgesehen.
# Aufgabe 3
Die Orientierung zu einer gegebenen Sequenz von Fragmenten kann berechnet werden, indem bei der Ausführung der Parameter **--consider_orientation** übergeben wird. Standardmäßig ist diese Option nun aktiviert, der Parameter muss also nicht zwingend übergeben werden.
Zur Funktionsweise: Zuerst wird die Orientierung berechnet und die daraus resultierende Orientierung wird als Input für den Sequenzer verwendet.
````bash
python assembler.py "data/frag.dat" --consider_orientation --verbose
````  
oder über die .exe:

````bash
assembler.exe "data/frag" --consider_orientation --verbose
````