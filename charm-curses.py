#!/usr/bin/env python

"""charm-curses.py: Simple command line interface for CHarm."""

import curses


class Screen:
    def __init__(self):
        self.stdscr = curses.initscr()
        self.stdscr.clear()
        curses.noecho()
        curses.cbreak()
        curses.curs_set(0)
        self.stdscr.keypad(1)
        curses.start_color()
        self.stdscr.refresh()

    def end_curses(self):
        curses.nocbreak()
        self.stdscr.keypad(0)
        curses.echo()
        curses.endwin()


class Window:
    def __init__(self, stdscr, height, width, x=0, y=0, box=False):
        self.width = width
        self.height = height
        self.x = x
        self.y = y
        self.box = box
        self.stdscr = stdscr

        self.window = curses.newwin(self.height, self.width, self.y, self.x)
        self.clean()


    def clean(self):
        self.window.clear()
        self.update()

    def print(self, string, x=1, y=1):
        self.window.addstr(y, x, string)
        self.update()

    def update(self):
        if self.box:
            self.window.box()
        self.window.refresh()
        self.stdscr.refresh()


class MessageBox(Window):
    def __init__(self, stdscr, box=True):
        self.width = width
        self.height = height
        self.x = x
        self.y = y
        self.box = box
        self.stdscr = stdscr

        self.window = curses.newwin(self.height, self.width, self.y, self.x)
        self.clean()


def main():
    screen = Screen()
    title_win = Window(screen.stdscr, 3, 80, 0, 0, box=True)
    main_win = Window(screen.stdscr, 21, 80, 0, 3, box=True)

    while True:
        screen.stdscr.clear()
        screen.stdscr.refresh()

        title_win.print("CHarm")
        main_win.update()

        c = screen.stdscr.getch()
        if c == ord('x'):
            break

    screen.end_curses()


if __name__ == "__main__":
    main()