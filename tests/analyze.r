--  This file is part of LS² - the Localization Simulation Engine of FU Berlin.
--
--  Copyright 2011-2013   Heiko Will, Marcel Kyas, Thomas Hillebrandt,
--  Stefan Adler, Malte Rohde, Jonathan Gunthermann
--
--  This program is free software: you can redistribute it and/or modify
--  it under the terms of the GNU General Public License as published by
--  the Free Software Foundation, either version 3 of the License, or
--  (at your option) any later version.
--
--  This program is distributed in the hope that it will be useful,
--  but WITHOUT ANY WARRANTY; without even the implied warranty of
--  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
--  GNU General Public License for more details.
--
--  You should have received a copy of the GNU General Public License
--  along with this program.  If not, see <http://www.gnu.org/licenses/>.
--

X <- read.csv("numbers.dat", head=FALSE)
X <- X[,1]

summary(X)
sd(X)
summary(abs(X))
sd(abs(X))
pdf("distribution.pdf", onefile=FALSE, height=4, width=6, pointsize=10)
hist(X[X<=500], prob=TRUE)
lines(density(X[X<=500], bw=0.1))

