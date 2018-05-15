# Install instructions (incomplete)

Install python3

Install rust at https://rustup.rs/

```
cd rust_fastforward
# This only ever needs to be run once, it installs build dependencies.
sudo pip3 install -r requirements.txt
# This needs to be run every time you change the rust code.
python3 setup.py install --user
```
