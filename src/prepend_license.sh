for f in *.cpp *.hpp *.c *.h; do
  cat ../license_header.txt $f > $f.new
  mv $f.new $f
done
