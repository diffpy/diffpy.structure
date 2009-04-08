PKGNAME = diffpy.Structure
OUTPUT = ../build/epydoc-html
WEBDOCDIR = login.cacr.caltech.edu:/cacr/home/proj/danse/docs.danse.us/docroot/diffraction


all:
	mkdir -p $(dir $(OUTPUT))
	epydoc --config=epydoc.cfg --output=$(OUTPUT)


clean:
	/bin/rm -rf $(OUTPUT)


publish:
	rsync -e /usr/bin/ssh -va --chmod=a+rwx,o-w,Dg+s,Fa-x \
	    $(OUTPUT)/ $(WEBDOCDIR)/$(PKGNAME)
