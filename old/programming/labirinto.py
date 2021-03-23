#!/usr/bin/env python
import numpy as np
import random

class Labirinto(object):
    """ An example maze class """
    def __init__(self, fixo=True, size=10):
        if fixo:
            size = 6
            self.M = LabA.copy()
        else:
            self.M = gera_labirinto(size=size)
        self.size = size
        # Where M has the structure
        #     0    1   2       3    4
        # M(LEFT, UP, RIGHT, DOWN, CHECK_IF_VISITED)
        # where LEFT means: wall at left, etc.

        # Mark everything as not visited
        self.M[:,:,4] = False
        # Except the very first position
        self.M[0,0,4] = True
        # Open the walls at the finish
        self.M[size-1,size-1,2] = True
        # Initializes the position and direction of the explorer
        self._x = 0
        self._y = 0
        self._direcao = 2
        self._livre_var = False

    def livre(self):
        """ Checks whether the explorer is free """
        return self._livre_var

    def movimenta(self, m):
        """ Moves the explorer.
            Input: m -> can be 'avançar' to move one cell
                   'direita' to turn around itself to the right
                   'esquerda' to turn around itself to the left
        """
        
        if m in ('avancar', 'avançar'):
            if self._livre_var:
                return True

            if self.M[self._y,self._x,self._direcao]:
                if self._direcao == 0:
                    self._x -= 1
                elif self._direcao == 1:
                    self._y -= 1
                elif self._direcao == 2:
                    if self._x == self.size-1 and self._y == self.size-1:
                        self._livre_var = True
                    self._x += 1
                elif self._direcao == 3:
                    self._y += 1

                if not self._livre_var:
                    self.M[self._y,self._x,4] = True

                return True
            else:
                return False
        elif m == 'direita':
            self._direcao +=1
            if self._direcao == 4:
                self._direcao = 0
            return True
        elif m == 'esquerda':
            self._direcao -=1
            if self._direcao == -1:
                self._direcao = 3
            return True
        else:
            raise ValueError

    def mostra(self, grava_no_arquivo=None):
        """ Shows the maze and the explorer """
        from matplotlib import pyplot as plt
        import matplotlib.cm as cm

        # Prepares the image containing the walls
        image = np.zeros((self.size*10,self.size*10), dtype=np.bool)
        for row in range(0,self.size):
            for col in range(0,self.size):
                cell_data = self.M[row,col]
                for i in range(10*row+1,10*row+9):
                    image[i,range(10*col+1,10*col+9)] = 255
                    if cell_data[0] == True:
                        image[range(10*row+1,10*row+9),10*col] = 255
                    if cell_data[1] == True:
                        image[10*row,range(10*col+1,10*col+9)] = 255
                    if cell_data[2] == True:
                        image[range(10*row+1,10*row+9),10*col+9] = 255
                    if cell_data[3] == True:
                        image[10*row+9,range(10*col+1,10*col+9)] = 255

        # Displays the walls
        ax = plt.subplot(111)
        plt.imshow(image, cmap = cm.Greys_r, interpolation='none',
                   extent=(0,self.size,self.size, 0))

        # Plots a red circle in visited cells
        for i in range(self.size):
            for j in range(self.size):
                if self._x == i and self._y == j:
                    continue
                if self.M[j,i,4]:
                    plt.plot(i+0.5, j+0.5, marker='o', color='r')

        # Plots the explorer, accounting for its direction
        marker_dict = { 2: '>', 0: '<', 1: '^', 3: 'v'}
        plt.plot(self._x+0.5, self._y+0.5,
                 marker=marker_dict[self._direcao],
                 color='b', markersize=14)

        # Sets grid and turns off tick labels
        plt.grid()
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        if grava_no_arquivo is not None:
            plt.savefig(grava_no_arquivo)
        plt.show()


def gera_labirinto(size=6):
    """ Generates a square maze using a random depth-first algorithm
        Input: optiona, size -> integer containing the side of the maze
                                Default: 6
        Returns: A NxNx5 numpy boolean array containing the details of the
                 maze. M[:,:,0] tells whether there is a wall

        Note: this function was adptaded from the maze-generation example
        algorithm in Wikipedia (retrieved Jul, 2016), originally coded by
        Erik Sweet and Bill Basener.
    """

    M = np.zeros((size,size,5), dtype=np.bool)
    # The array M is going to hold the array information for each cell.
    # The first four coordinates tell if walls exist on those sides
    # and the fifth indicates if the cell has been visited in the search.
    # M(LEFT, UP, RIGHT, DOWN, CHECK_IF_VISITED)
    image = np.zeros((size*10,size*10), dtype=np.bool)
    # The array image is going to be the output image to display

    # Set starting row and column
    r = 0
    c = 0
    history = [(r,c)] # The history is the

    # Trace a path though the cells of the maze and open walls along the path.
    # We do this with a while loop, repeating the loop until there is no history,
    # which would mean we backtracked to the initial start.
    while history:
        M[r,c,4] = True # designate this location as visited
        # check if the adjacent cells are valid for moving to
        check = []
        if c > 0 and M[r,c-1,4] == False:
            check.append('L')
        if r > 0 and M[r-1,c,4] == False:
            check.append('U')
        if c < size-1 and M[r,c+1,4] == False:
            check.append('R')
        if r < size-1 and M[r+1,c,4] == False:
            check.append('D')

        if len(check): # If there is a valid cell to move to.
            # Mark the walls between cells as open if we move
            history.append([r,c])
            move_direction = random.choice(check)
            if move_direction == 'L':
                M[r,c,0] = True
                c = c-1
                M[r,c,2] = True
            if move_direction == 'U':
                M[r,c,1] = True
                r = r-1
                M[r,c,3] = True
            if move_direction == 'R':
                M[r,c,2] = True
                c = c+1
                M[r,c,0] = True
            if move_direction == 'D':
                M[r,c,3] = True
                r = r+1
                M[r,c,1] = True
        else: # If there are no valid cells to move to.
        # retrace one step back in history if no move is possible
            r,c = history.pop()
    return M


# An example maze (generated with gera_labirinto and selected by hand)
LabA = np.array([[[False, False,  True, False,  True,],
  [ True, False, False,  True,  True,],
  [False, False,  True,  True,  True,],
  [ True, False,  True, False,  True,],
  [ True, False,  True, False,  True,],
  [ True, False, False,  True,  True,],],

 [[False, False,  True,  True,  True,],
  [ True,  True, False, False,  True,],
  [False,  True, False,  True,  True,],
  [False, False,  True,  True,  True,],
  [ True, False, False, False,  True,],
  [False,  True, False,  True,  True,],],

 [[False,  True,  True, False,  True,],
  [ True, False, False,  True,  True,],
  [False,  True, False,  True,  True,],
  [False,  True,  True,  True,  True,],
  [ True, False,  True, False,  True,],
  [ True,  True, False, False,  True,],],

 [[False, False, False,  True,  True,],
  [False,  True, False,  True,  True,],
  [False,  True, False,  True,  True,],
  [False,  True,  True, False,  True,],
  [ True, False, False,  True,  True,],
  [False, False, False,  True,  True,],],

 [[False,  True, False,  True,  True,],
  [False,  True,  True, False,  True,],
  [ True,  True, False, False,  True,],
  [False, False,  True,  True,  True,],
  [ True,  True, False, False,  True,],
  [False,  True, False,  True,  True,],],

 [[False,  True,  True, False,  True,],
  [ True, False,  True, False,  True,],
  [ True, False,  True, False,  True,],
  [ True,  True,  True, False,  True,],
  [ True, False,  True, False,  True,],
  [ True,  True, False, False,  True,],],])
