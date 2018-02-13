"""
IGV related code.
"""
import os
from recipes.render import render_template
from django import template

register = template.Library()

# Known filetypes.
KNOWN_EXT = ["bam", "bw", "vcf",
             "gff3", "gtf", "bedgraph", "bigwig"]

class Bunch(object):
    def __init__(self, path, baseurl=''):
        self.path = path
        self.name = os.path.basename(self.path)
        self.link = os.path.join(baseurl, self.name)

@register.inclusion_tag('igv/bam.xml')
def bam(bunch):
    """
    Generates tracks for bam files.
    """
    return dict(bunch=bunch)


@register.inclusion_tag('igv/bigwig.xml')
def bigwig(path, name):
    """
    Generates tracks for coverage files.
    """
    return dict(path=path, name=name)


@register.inclusion_tag('igv/resources.xml')
def resources(path):
    """
    Generates a resource tag.
    """
    return dict(path=path)


def walk(root, baseurl, collect=[]):
    for entry in os.scandir(root):
        if entry.is_file():
            bunch = Bunch(path=entry.path, baseurl=baseurl)
            collect.append(bunch)
        else:
            walk(root=entry.path, baseurl=baseurl, collect=collect)
    return collect


def create_resources(root, baseurl):
    data = list(walk(root, baseurl=baseurl))
    return data

def filter_resources(data, type='bam'):
    subset = filter(lambda bunch: bunch.path.endswith(type), data)
    subset = list(subset)
    return subset

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(description='Creates an IGV session '
                                                 'file from all files in a directory.')

    parser.add_argument('--root', default="", help='The location that needs to be searched')
    parser.add_argument('--url', help='The base URL for the directory')
    parser.add_argument('--genome', default="genome.fa", help='The path to the genome')

    return parser


if __name__ == '__main__':

    genome="/Users/ialbert/edu/24/db/AF086833.fa"

    data = create_resources(root="/Users/ialbert/edu/24/", baseurl='/Users/ialbert/edu/24/')
    bams = filter_resources(data, 'bam')

    context = dict(bams=bams, genome=genome)

    html = render_template("igv/main.xml", context=context)

    print (html)
