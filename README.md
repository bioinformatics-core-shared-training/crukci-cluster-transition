# CRUKCI Cluster Transition - Hands-on training sessions

4 hands-on training sessions for using the cluster to run analysis @ CRUKCI for complete novice.

## What do you need?

- A **computer** (ideally a Unix-like one) connected to the Internet
- A **cluster account** -- request one via Helpdesk - IT <ithelpdesk@cruk.cam.ac.uk>
- A **text editor** -- download [Atom](https://atom.io/) if you do not have one already
- A **web browser** to access [this page](https://github.com/bioinformatics-core-shared-training/crukci-cluster-transition)

and the ability to **submit jobs onto the cluster**, please check if you can by following instructions on [Can I submit jobs onto the cluster?](can-i-submit-jobs.md) before doing sessions 3 & 4.

## Where to start?

- [Session 1: Shell](session1-shell.md)
- [Session 2: Cluster](session2-cluster.md)
- [Session 3: Analysis steps by steps](session3-analysis.md)
- [Session 4: Advanced analysis and usage of the cluster](session4-advanced.md)

## Tiny url of this page

- [https://frama.link/crukci-cluster](https://frama.link/crukci-cluster)

## Rendering GitHub markdown files locally

https://github.com/joeyespo/grip

```shell
python3 -m venv venv
source venv/bin/activate
pip install grip
grip session1-shell.md
```
