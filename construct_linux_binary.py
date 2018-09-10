"""
This script copies what the Travis build is doing, but allows you to build in your own environment
"""

import argparse
import subprocess
import os
import sys
import github3
import shutil

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument('--token', required=True)
    args.add_argument('--repo', default='vartrix')
    args.add_argument('--owner', default='10xgenomics')
    args.add_argument('--target', default='target/release/vartrix')
    args.add_argument('--extra-items', nargs='+', default=['README.md', 'LICENSE'])
    args.add_argument('--tag', help='Will use the last release if not set')
    return args.parse_args()


def construct_tarball(target, extra_items, file_name):
    """Package a tarball of the binary + any extra items"""
    o, e = subprocess.Popen(['mktemp'], stdout=subprocess.PIPE).communicate()
    tar_dir = file_name.replace('.tar.gz', '')
    os.makedirs(os.path.join(o, tar_dir))
    if e:
        raise Exception('failed to make temp dir')
    shutil.copy(target, os.path.join(o, tar_dir, os.path.basename(target)))
    for f in extra_items:
        shutil.copy(f, os.path.join(o, tar_dir, os.path.basename(f)))
    subprocess.check_call(['strip', os.path.join(o, tar_dir, os.path.basename(target))])
    os.chdir(o)
    subprocess.check_call(['tar', 'czf', file_name, tar_dir])
    tarball = os.path.join(o, file_name)
    os.chdir('../')
    return tarball


def upload_to_github(args):
    github = github3.login(token=args.token)
    repo = github.repository(args.owner, args.repo)
    releases = list(repo.releases())
    
    if args.tag is not None:
        tag = args.tag
        try:
            release = [x for x in releases if x.tag_name == args.tag][0]
        except IndexError:
            print "WARNING: provided tag does not exist. Using last release"
            release = releases[-1]
    else:
        release = releases[-1]
        tag = release.tag_name
    
    file_name = 'vartrix-{}-x86_64-linux.tar.gz'.format(tag)
    tarball = construct_tarball(args.target, args.extra_items, file_name)
    try:
        release.upload_asset('application', file_name, open(tarball, 'rb'), label=file_name)
    except github3.GitHubError as e:
        raise Exception(e.errors)


if __name__ == '__main__':
    args = parse_args()
    upload_to_github(args)