# Foreword

The extrafont R package was created and maintained by [Winston Chang](https://github.com/wch) on [github](https://github.com/wch/extrafont) for more than 13 years with a first commit on github the 17th of May 2012.
The previous commit history of the package can be found [here](https://github.com/wch/extrafont).
The responsibilities and attention of Winston Chang are on things quite far away from extrafont and Rttf2pt1 these days. 
As a consequence he asked for another person to maintain the extrafont package and I ([Frederic Bertrand](https://github.com/fbertran)) volunteered [(read here)](https://github.com/wch/Rttf2pt1/issues/25#issuecomment-3320579566).

The following text is the original readme from Winston Chang repository on [github](https://github.com/wch/extrafont), only some links were updated.

-------------------

# extrafont

The extrafont package makes it easier to use fonts other than the basic PostScript fonts that R uses.
Fonts that are imported into extrafont can be used with PDF or PostScript output files. On Windows, extrafont will also make system fonts available for bitmap output.

There are two hurdles for using fonts in PDF (or Postscript) output files:

* Making R aware of the font and the dimensions of the characters.
* Embedding the fonts in the PDF file so that the PDF can be displayed properly on a device that doesn't have the font. This is usually needed if you want to print the PDF file or share it with others.

The extrafont package makes both of these things easier.

Presently it allows the use of TrueType fonts with R, and installation of special font packages.
Support for other kinds of fonts will be added in the future.
It has been tested on Mac OS X 10.7 and Ubuntu Linux 12.04 and Windows XP.


The instructions below are written for PDF files, although the information also applies to PostScript files.

If you want to use the TeX Computer Modern fonts in PDF files, also see the [fontcm](https://github.com/wch/fontcm) package.

# Using extrafont

## Requirements

You must have Ghostscript installed on your system for embedding fonts into PDF files.

Extrafont requires the **[extrafontdb](https://github.com/wch/extrafontdb)** package to be installed.
extrafontdb contains the font database, while this package contains the code to install fonts and register them in the database.

It also requires the **[Rttf2pt1](https://github.com/fbertran/Rttf2pt1)** package to be installed.
Rttf2pt1 contains the ttf2pt1 program which is used to read and manipulate TrueType fonts.
It is in a separate package for licensing reasons.

Install extrafont from CRAN will automatically install extrafontdb and Rttf2pt1:

```R
install.packages('extrafont')
library(extrafont)
```


To use extrafont in making graphs, you'll need to do the following:

* Import fonts into the extrafont database. (Needs to be done once)
* Register the fonts from the extrafont database with R's PDF (or PostScript) output device. (Needs to be done once per R session)
* Create the graphics that use the fonts.
* Embed the fonts into the PDF file. (Needs to be done for each file)

## Import fonts into the extrafont database

First, import the fonts installed on the system.
(This only works with TrueType fonts right now.)

```R
font_import()
# This tries to autodetect the directory containing the TrueType fonts.
# If it fails on your system, please let me know.
```

This does the following:

* Finds the fonts on your system.
* Extracts the FontName (like ArialNarrow-BoldItalic).
* Extracts/converts a PostScript .afm file for each font. This file contains the *font metrics*, which are the rectangular dimensions of each character that are needed for placement of the characters. These are not the *glyphs*, which the curves defining the visual shape of each character. The glyphs are only in the .ttf file.
* Scan all the resulting .afm files, and save a table with information about them. This table will be used when making plots with R.
* Creates a file `Fontmap`, which contains the mapping from FontName to the .ttf file. This is required by Ghostscript for embedding fonts.


You can view the resulting table of font information with:

```R
# Vector of font family names
fonts()

# Show entire table
fonttable()
```

If you install new fonts on your computer, you'll have to run `font_import()` again.

## Register the fonts with the PDF output device

The next step is to register the fonts in the afm table with R's PDF (or PostScript) output device.
This is needed to create PDF files with the fonts.
As of extrafont version 0.13, this must be run only in the first session when you import your fonts.
In sessions started after the fonts have been imported, simply loading the package with `library(extrafont)` this step isn't necessary, since it will automatically register the fonts with R.

```R
# Only necessary in session where you ran font_import()
loadfonts()
# For PostScript output, use loadfonts(device="postscript")
# Suppress output with loadfonts(quiet=TRUE)
```


## Create figures with the fonts


Here's an example of PDFs made with base graphics and with ggplot2.
These examples use the font Impact, which should be available on Windows and Mac.
(Use `fonts()` to see what fonts are available on your system)

```R
pdf("font_plot.pdf", family="Impact", width=4, height=4)
plot(mtcars$mpg, mtcars$wt, 
     main = "Fuel Efficiency of 32 Cars",
     xlab = "Weight (x1000 lb)",
     ylab = "Miles per Gallon")
dev.off()


library(ggplot2)
p <- ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point() +
  ggtitle("Fuel Efficiency of 32 Cars") +
  xlab("Weight (x1000 lb)") + ylab("Miles per Gallon") +
  theme(text=element_text(size=16, family="Impact"))

ggsave("font_ggplot.pdf", plot=p,  width=4, height=4)
```

The first time you use a font, it may throw some warnings about unknown characters.
This should be harmless, but if it causes any problems, please report them.


## Embed the fonts

After you create a PDF output file, you should *embed* the fonts into the file.
There are 14 PostScript *base fonts* never need to be embedded, because they are included with every PDF/PostScript renderer.
All other fonts should be embedded into the PDF files.


First, if you are running Windows, you may need to tell it where the Ghostscript program is, for embedding fonts. (See Windows installation notes below.)

```R
# Needed only on Windows - run once per R session
# Adjust the path to match your installation of Ghostscript
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.05/bin/gswin32c.exe")
```


As the name suggests, `embed_fonts()` will embed the fonts:

```R
embed_fonts("font_plot.pdf", outfile="font_plot_embed.pdf")
embed_fonts("font_ggplot.pdf", outfile="font_ggplot_embed.pdf")
# If outfile is not specified, it will overwrite the original file
```

To check if the fonts have been properly embedded, open each of the PDF files with Adobe Reader, and go to File->Properties->Fonts.
If a font is embedded, it will say "Embedded Subset" by the font's name; otherwise it will say nothing next to the name.

With Adobe Reader, if a font is not embedded, it will be substituted by another font.
This provides a way to see what your PDF will look like on printer or computer that doesn't have the font installed.
Other PDF viewers may behave differently.
For example, the Preview application on Mac OS X will automatically use system fonts to display non-embedded fonts -- this makes it impossible to tell whether the font is embedded in the PDF.

On Linux you can also use evince (the default PDF viewer) to view embedded fonts.
Open the file and go to File->Properties->Fonts.
If a font is embedded, it will say "Embedded subset"; otherwise it will say "Not embedded".

If you are putting multiple PDF figures into a single document, it is more space-efficient to _not_ embed fonts in each figure, but instead embed the font in the final PDF document.

## Windows bitmap output

extrafont also makes it easier to use fonts in Windows for on-screen or bitmap output.

```R
# Register fonts for Windows bitmap output
loadfonts(device="win")

ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point() +
  ggtitle("Title text goes here") +
  theme(plot.title = element_text(size = 16, family="Georgia", face="italic"))

ggsave("fonttest-win.png")
```

Since the output is a bitmap file, there's no need to embed the fonts.

# Font packages

Extrafont supports _font packages_, which contain fonts that are packaged in a particular way so that they can be imported into extrafont.
These fonts are installed as R packages; they are not installed for the computer operating system.
Fonts that are installed this way will be available only for PDF or PostScript output.
They will not be available for on-screen or bitmap output, which requires that the font be installed for operating system, not just with R and extrafont.

Presently extrafont supports only font packages with PostScript Type 1 fonts.

See the [fontcm](https://github.com/wch/fontcm) package containing Computer Modern fonts for an example.

*****

# Installation notes

## Rttf2pt1

The source code for the utility program `ttf2pt1` is in the package Rttf2pt1.
CRAN has pre-compiled Windows and Mac OS X binaries.
For other platforms, and when installing from source, it will be compiled on installation, so you need a build environment on your system.


## Windows installation notes

In Windows, you need to make sure that Ghostscript is installed.

In each R session where you embed fonts, you will need to tell R where Ghostscript is installed.
For example, when Ghostscript 9.05 is installed to the default location, running this command will do it (adjust the path for your installation):

```R
Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.05/bin/gswin32c.exe")
```

## Resetting the font database

To reset the extrafont database, reinstall the extrafontdb package:

```R
install.packages("extrafontdb")
```
