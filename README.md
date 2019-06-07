# TM25RaySetTools
A set of C++ tools for IES-TM25 ray sets
Created by Julius Muschaweck, JMO GmbH, Gauting, Germany (julius@jmoptics.de)

Background: In illumination optics, ray files are often used as source models
in Monte Carlo simulations.
Ray files are typically large, containing millions of random rays. Traditionally,
each optical simulation software vendor uses his proprietary file format, which
is not only a web space and quality control nightmare for the LED vendors. The 
proprietary formats are also often poorly documented and not flexible.
Therefore, an IES committee created a standardized, well engineered, flexible,
compact and well documented file format, containing a superset of all features of all 
major proprietary formats. This file format is called IES-TM25-13, or TM25 for short.
The standard, published in 2013, is available at
https://www.ies.org/product/ray-file-format-for-the-description-of-the-emission-property-of-light-sources/ 

TM25 ray files are now becoming more common. Simulation software vendors have
implemented the ability to read TM25 ray files, and LED vendors provide them.
For my own scientific and engineering work, I need direct access to TM25 ray file
content. Therefore, I wrote C++ software to parse TM25 ray files and provide the
content in a C++ class, and I make this software available to the community under 
the "unlicense", see unlicense.txt, which basically allows you to do whatever you 
want with the code, except sueing me if it doesn't work for you.

Currently, the software is available as a Visual Studio 2017 project. I have tried
my best to write it in portable way, using only C++ constructs from the C++11 standard.
I'm planning to test it on a Linux/gcc platform as well, and to add more features.

As always, I'm looking forward to your feedback. Please use the standard GitHub workflow,
adding issues, and pull requests if you'd like to provide your own additions / bug fixes.

Enjoy!

Gauting, Germany, on a sunny day in 2019
Julius Muschaweck
