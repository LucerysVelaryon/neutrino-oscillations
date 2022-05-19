import time
import sys

class BarreDeProgression:
    """Cette classe permet de créer des barres de chargement"""

    def __init__(self, taille=30, titre='Chargement'):
        self.taille = taille
        self.titre = titre
        self.pourcentage = 0
        self.pourcentage_precedent = 0
        self.temps = time.time()
        self.temps_restant = 'Calcul...'

        print("\nC'est parti !")
        self.maj()

    def maj(self, pourcentage=0):
        self.pourcentage = pourcentage
        etapes = int(self.pourcentage / 100 * self.taille)

        if self.pourcentage - self.pourcentage_precedent >= 5:
            self.temps_restant = (time.time() - self.temps) * (100 - pourcentage) / 5
            self.temps = time.time()
            self.pourcentage_precedent = self.pourcentage

        temps_restant = self.temps_restant

        if etapes == 0:
            visuel = self.taille * ' '
        else:
            visuel = etapes * '=' + (self.taille - etapes) * ' '

        if self.pourcentage == 100:
            sys.stdout.write('\rTerminé !' + (self.taille + 100) * ' ' + '\n')
        else:
            if type(self.temps_restant) != str:
                temps_restant = format_temps(int(self.temps_restant))
            sys.stdout.write('\r' + self.titre + ' [' + visuel + '] ' + str(
                int(self.pourcentage)) + '%        Temps restant : ' + temps_restant)
        sys.stdout.flush()


def format_temps(temps):
    if temps > 3600:
        h = temps // 3600
        m = (temps // 60) % 60
        s = temps % 60
        return str(h) + ' h. ' + str(m) + ' m. ' + str(s) + ' s.   '
    elif temps > 60:
        m = temps // 60
        s = temps % 60
        return str(m) + ' m. ' + str(s) + ' s.          '
    else:
        return str(temps) + ' s.                '
